using System.Collections.Concurrent;

namespace BinChillingTools;


public sealed class BinRefiner
{
    private readonly int _partitionSizeUpperBound;
    private readonly IReadOnlyDictionary<CoTuple, double> _coMatrix;
    private readonly BinEvaluator _evaluator;
    private readonly ConcurrentDictionary<object, double> _scoreDict = new();


    public BinRefiner(int partitionSizeUpperBound, IReadOnlyDictionary<CoTuple, double> coMatrix, BinEvaluator evaluator)
    {
        _partitionSizeUpperBound = partitionSizeUpperBound;
        _coMatrix = coMatrix;
        _evaluator = evaluator;
    }

    private readonly record struct CoResult(double Score, string[] Cluster);


    public IReadOnlyList<string[]> Refine(IReadOnlyList<string[]> partition)
    {
        Console.WriteLine($"Starting refinement using {partition.Count} clusters");
        
        _scoreDict.Clear();
        var initialClusters = partition as string[][] ?? partition.ToArray();
        
        var removeSet = new HashSet<object>();
        
        double internalCoScore;
        using (var bar = ProgressBarHandler.CreateBar(initialClusters.Length))
        {
            internalCoScore = initialClusters.InternalCoScoreParallel(_coMatrix, bar);
            bar.Final();
        }
        Console.WriteLine($"Average co score: {internalCoScore}");

        var currentSize = initialClusters.Length;
        var currentSegment = initialClusters.Segment();

        var oldClusters = new List<string[]>();
        List<string[]> newClusters;
        
        while (true)
        {
            using (var progressBar = ProgressBarHandler.CreateBar(currentSize))
            {
                newClusters = RefineClustersMultiThreaded_NoBreak_NoRepeat(
                    currentSegment.UseProgressBar(progressBar), oldClusters, internalCoScore, removeSet);
                progressBar.Final();
            }
            
            
            oldClusters.RemoveAll(removeSet.Contains);
            Console.WriteLine($"Refined {removeSet.Count} clusters into {newClusters.Count} new clusters." +
                              $" ({oldClusters.Count + newClusters.Count})");
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
        var result = newClusters.Concat(oldClusters).ToList();
        Console.WriteLine($"Average co score: {result.InternalCoScoreParallel(_coMatrix)}");
        return result;
    }


    private double GetClusterScoreDynamic(IEnumerable<string> cluster)
    {
        if(_scoreDict.TryGetValue(cluster, out var score)) return score;
        score = _evaluator.ScoreCluster(cluster);
        _scoreDict[cluster] = score;
        return score;
    }
    

    private List<string[]> RefineClustersMultiThreaded_NoBreak_NoRepeat(
        IEnumerable<KeyValuePair<string[], IReadOnlyCollection<string[]>>> clusterList, ICollection<string[]> oldClusters,
        double internalCoScore, ISet<object> skipSet)
    {
        var newClusters = new List<string[]>();
        foreach (var (cluster1, clusterList2) in clusterList)
        {
            if(skipSet.Contains(cluster1) || cluster1.Length == 0) continue;
            var score1 = GetClusterScoreDynamic(cluster1);

            CoResult? scoreResult = null;
            scoreResult = SearchOptimalPair(cluster1, clusterList2, score1, internalCoScore, skipSet);

            //Add best cluster, if any
            if (scoreResult.HasValue)
            {
                skipSet.Add(cluster1);
                skipSet.Add(scoreResult.Value.Cluster);
                newClusters.Add( cluster1.Concat(scoreResult.Value.Cluster).ToArray() );
                continue;
            }

            if (score1 < 0.0d)
            {
                skipSet.Add(cluster1);
                newClusters.AddRange( SplitBin( cluster1 ) );
                continue;
            }
            
            //only added if not merged with another cluster
            oldClusters.Add(cluster1);
        }
        return newClusters;
    }

    private CoResult? SearchOptimalPair(IReadOnlyList<string> cluster1, IEnumerable<string[]> clusterList2,
        double score1, double internalCoScore, ISet<object> skipSet)
    {
        CoResult? scoreResult = null;
        Parallel.ForEach(clusterList2
                .Where(c => !skipSet.Contains(c) && c.Length > 0),
            cluster2 =>
            {
                //only calc score2 if score1 exists.
                var score2 = double.IsNaN(score1)
                    ? double.NaN
                    : GetClusterScoreDynamic(cluster2);
                
                double delta, totalScore;
                //use scg score if both score1 and score2 has it.
                if (double.IsNaN(score1) || double.IsNaN(score2))
                {
                    var commonCo = CommonCoCalculation(cluster1, cluster2, _coMatrix);
                    delta = commonCo - internalCoScore;
                    if (delta <= 0) return;
                    totalScore = commonCo;
                }
                else
                {
                    //calculate score first, as is typically faster and might fail more often
                    if (double.IsNaN(score2)) return;
                    var comboScore = _evaluator.ScoreCluster(cluster1.Concat(cluster2));
                    delta = comboScore - Math.Max(score1, score2);
                    if (delta <= 0) return;

                    var commonCo = CommonCoCalculation(cluster1, cluster2, _coMatrix);
                    totalScore = commonCo * (1+comboScore);
                }

                //not locked, as to have maximum throughput
                if (scoreResult is null || totalScore > scoreResult.Value.Score)
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
        for (var i = 0; i < cluster1.Count; i++)
        {
            for (var j = 0; j < cluster2.Count; j++)
            {
                var val = coMatrix.TryGetValue(new CoTuple(cluster1[i], cluster2[j]), out var value ) ? value : 0.0d;
                sum = sum + val;
            }
        }

        return sum / (cluster1.Count * cluster2.Count);
    }
    
}

