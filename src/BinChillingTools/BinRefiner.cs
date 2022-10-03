using System.Collections.Concurrent;
using ShellProgressBar;

namespace BinChillingTools;


public sealed class BinRefiner
{
    private readonly int _partitionSizeUpperBound;
    private readonly double _minCoValue;
    private readonly IReadOnlyDictionary<CoTuple, double> _coMatrix;
    private readonly BinEvaluator _evaluator;
    private readonly ConcurrentDictionary<int, double> _scoreDict = new();

    public BinRefiner(int partitionSizeUpperBound, double minCoValue, IReadOnlyDictionary<CoTuple, double> coMatrix, BinEvaluator evaluator)
    {
        _partitionSizeUpperBound = partitionSizeUpperBound;
        _minCoValue = minCoValue;
        _coMatrix = coMatrix;
        _evaluator = evaluator;
    }

    private readonly record struct RefineResult(List<string[]> OldClusters, List<string[]> NewClusters);
    private readonly record struct CoResult(double Score, string[] Cluster);



    public IReadOnlyList<string[]> Refine(IReadOnlyDictionary<string, string[]> partition)
    {
        Console.WriteLine($"Starting refinement using {partition.Count} clusters");
        
        _scoreDict.Clear();
        var namedClusters = partition
            .Select(x => x.Value)
            .ToArray();
        
        var removeSet = new HashSet<int>();

        var currentSize = namedClusters.Length;
        var currentSegment = namedClusters.Segment();
        var oldClusters = new List<string[]>();
        List<string[]> newClusters;


        while (true)
        {
            newClusters = RefineClustersMultiThreaded_NoBreak_NoRepeat(currentSegment, oldClusters, currentSize, removeSet);
            Console.WriteLine($"Refined {removeSet.Count} clusters into {newClusters.Count} new clusters.");

            
            oldClusters.RemoveAll(x => removeSet.Contains(x.GetHashCode()));
            foreach (var clusterHash in removeSet) 
                _scoreDict.Remove(clusterHash, out _);

            if (newClusters.Count < 2 || (oldClusters.Count + newClusters.Count) <= _partitionSizeUpperBound ) 
                break;
            
            //update current iterator
            currentSize = newClusters.Count;
            currentSegment = newClusters.SmartGroupSegment(oldClusters);
            
            //reset removeSet
            removeSet.Clear();
            GC.Collect();
        }
        
        _scoreDict.Clear();
        return newClusters.Concat(oldClusters).ToList();
    }


    private double GetClusterScoreDynamic(IEnumerable<string> cluster)
    {
        var hash = cluster.GetHashCode();
        if(_scoreDict.TryGetValue(hash, out var score)) return score;
        score = _evaluator.ScoreCluster(cluster);
        // Console.WriteLine($"score: {score}");
        _scoreDict[hash] = score;
        return score;
    }

    private List<string[]> RefineClustersMultiThreaded_NoBreak_NoRepeat(
        IEnumerable<KeyValuePair<string[], IReadOnlyCollection<string[]>>> clusterList, ICollection<string[]> oldClusters, int size, ISet<int> skipSet)
    {
        using var progressBar = ProgressBarHandler.Create(size);
        var newClusters = new List<string[]>();

        foreach (var (cluster1, clusterList2) in clusterList.UseProgressBar(progressBar))
        {
            if(skipSet.Contains(cluster1.GetHashCode()) || cluster1.Length == 0) continue;
            var score1 = GetClusterScoreDynamic(cluster1);

            CoResult? scoreResult = null;
            using (var childBar = progressBar.Spawn(clusterList2.Count, ProgressBarHandler.DefaultStartMsg,
                       ProgressBarHandler.DefaultChildOptions))
            {
                scoreResult = SearchOptimalPair(cluster1, clusterList2, score1, skipSet, childBar);
                childBar.Finalize();
            }

            //Add best cluster, if any
            if (scoreResult.HasValue)
            {
                skipSet.Add(cluster1.GetHashCode());
                skipSet.Add(scoreResult.Value.Cluster.GetHashCode());
                newClusters.Add( cluster1.Concat(scoreResult.Value.Cluster).ToArray() );
                continue;
            }

            if (score1 < 0.0d)
            {
                skipSet.Add(cluster1.GetHashCode());
                newClusters.AddRange( SplitBin( cluster1 ) );
                continue;
            }
            
            //only added if not merged with another cluster
            oldClusters.Add(cluster1);
        }
        progressBar.Finalize();
        return newClusters;
    }

    private CoResult? SearchOptimalPair(IReadOnlyList<string> cluster1, IReadOnlyCollection<string[]> clusterList2,
        double score1, ISet<int> skipSet, IProgressBar? progressBar = null)
    {
        var enumerable = progressBar is null ? clusterList2 : clusterList2.UseProgressBar(progressBar);
        CoResult? scoreResult = null;
        Parallel.ForEach(enumerable
                .Where(c => !skipSet.Contains(c.GetHashCode()) && c.Length > 0),
            cluster2 =>
            {
                var commonCo = CommonCoCalculation(cluster1, cluster2, _coMatrix);
                
                if (commonCo <= _minCoValue) return;
                var score2 = GetClusterScoreDynamic(cluster2);

                var comboScore = _evaluator.ScoreCluster(cluster1.Concat(cluster2));
                var delta = comboScore - Math.Max(score1, score2);
                var totalScore = comboScore * commonCo;
                
                //not locked, as to have maximum throughput
                if (delta > 0 && (scoreResult is null || totalScore > scoreResult.Value.Score))
                {
                    //lock skipset to go in here. Cant lco scoreResult, so using skipSet instead.
                    lock (skipSet)
                    {
                        //repeated as to not override good result
                        if (scoreResult is null || totalScore > scoreResult.Value.Score)
                        {
                            scoreResult = new CoResult(totalScore, cluster2);
                        }
                    }
                }
            });
        
        return scoreResult;
    }

    private static IEnumerable<string[]> SplitBin(IEnumerable<string> cluster)
    {
        return cluster.Select(x => new []{ x });
    }

    private static double CommonCoCalculation(IReadOnlyList<string> cluster1, IReadOnlyList<string> cluster2,
        IReadOnlyDictionary<CoTuple, double> coMatrix)
    {
        var sum = 0.0d;
        if (cluster1.Count == 0 || cluster2.Count == 0) return sum;
        for (int i = 0; i < cluster1.Count; i++)
        {
            for (int j = 0; j < cluster2.Count; j++)
            {
                var val = coMatrix.TryGetValue(new CoTuple(cluster1[i], cluster2[j]), out var value ) ? value : 0.0d;
                sum = sum + val;
            }
        }

        return sum / (cluster1.Count * cluster2.Count);
    }
    
}

