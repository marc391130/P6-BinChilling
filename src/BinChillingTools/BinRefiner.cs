using System.Collections.Concurrent;
using System.Diagnostics.CodeAnalysis;

namespace BinChillingTools;


public sealed class BinRefiner
{
    private readonly int _partitionSizeUpperBound;
    private readonly IReadOnlyDictionary<CoTuple, double> _coMatrix;
    private readonly BinEvaluator _evaluator;
    private readonly ConcurrentDictionary<object, double> _scoreDict = new();
    private readonly ConcurrentDictionary<object, double> _icoScoreDict = new();


    public BinRefiner(int partitionSizeUpperBound, IReadOnlyDictionary<CoTuple, double> coMatrix, BinEvaluator evaluator)
    {
        _partitionSizeUpperBound = partitionSizeUpperBound;
        _coMatrix = coMatrix;
        _evaluator = evaluator;
    }

    private readonly record struct CoResult(double Score, string[] Cluster);


    [SuppressMessage("ReSharper", "PossibleMultipleEnumeration")]
    public IReadOnlyList<string[]> Refine(IReadOnlyList<string[]> partition)
    {
        Console.WriteLine($"Starting refinement using {partition.Count} clusters");
        
        _scoreDict.Clear();
        // var initialClusters = SplitAllNegativeClusters(partition);
        var initialClusters = partition as List<string[]> ?? partition.ToList();
        var removeSet = new HashSet<object>();
        
        var oldClusters = new List<string[]>();
        var newClusters = initialClusters;
        initialClusters = null; //kill, so theres more memory
        GC.Collect();

        while (true)
        {
            using (var progressBar = ProgressBarHandler.CreateBar(newClusters.Count))
            {
                var segment = newClusters.Segment().ExtendSegment(oldClusters).UseProgressBar(progressBar);
                newClusters = RefineClustersMultiThreaded_NoBreak_NoRepeat(segment, oldClusters, removeSet);
                progressBar.Final();
            }
            
            
            oldClusters.RemoveAll(removeSet.Contains);
            // Console.WriteLine($"ICO: {internalCoScore}");
            Console.WriteLine($"Refined {removeSet.Count} clusters into {newClusters.Count} new clusters." +
                              $" (total: {oldClusters.Count + newClusters.Count})");
            foreach (var bin in removeSet)
            {
                _scoreDict.Remove(bin, out _);
                _icoScoreDict.Remove(bin, out _);

            }

            if (newClusters.Count == 0) break;

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
    
    private double GetIcoScoreDynamic(IReadOnlyList<string> cluster)
    {
        if(_icoScoreDict.TryGetValue(cluster, out var score)) return score;
        score = cluster.InternalCo(_coMatrix);
        _icoScoreDict[cluster] = score;
        return score;
    }
    

    private List<string[]> RefineClustersMultiThreaded_NoBreak_NoRepeat(
        IEnumerable<KeyValuePair<string[], IReadOnlyList<string[]>>> clusterList, ICollection<string[]> oldClusters, 
        ISet<object> skipSet)
    {
        var newClusters = new List<string[]>();
        foreach (var (cluster1, clusterList2) in clusterList)
        {
            if(skipSet.Contains(cluster1) || cluster1.Length == 0) continue;

            var score1IcoTask = Task.Factory.StartNew(() => GetIcoScoreDynamic(cluster1));
            var score1Task = Task.Factory.StartNew(() => GetClusterScoreDynamic(cluster1));

            var score1 = score1Task.Result; 
            var score1Ico = score1IcoTask.Result;
            
            var results = SearchOptimalPair(cluster1, clusterList2, score1, score1Ico, skipSet);

            //Add best cluster, if any
            if (results.HasValue)
            {
                skipSet.Add(cluster1);
                skipSet.Add(results.Value.Cluster);
                newClusters.Add( cluster1.Concat(results.Value.Cluster).ToArray() );
                continue;
            }

            if (!double.IsNaN(score1) && score1 < 0)
            {
                skipSet.Add(cluster1);
                newClusters.AddRange(SplitBin(cluster1));
                continue;
            }
            
            //only added if not merged with another cluster
            oldClusters.Add(cluster1);
        }
        return newClusters;
    }

    private CoResult? SearchOptimalPair(IReadOnlyList<string> cluster1, IEnumerable<string[]> clusterList2,
        double score1, double score1Ico, ISet<object> skipSet)
    {
        var useScGs = false;
        CoResult? scoreResult = null;
        Parallel.ForEach(clusterList2
                .Where(c => !skipSet.Contains(c) && c.Length > 0),
            cluster2 =>
            {
                //only calc score2 if score1 exists.
                var localScgUsage = false;
                var score2 = double.IsNaN(score1)
                    ? GetIcoScoreDynamic(cluster2)
                    : GetClusterScoreDynamic(cluster2);
                
                //use scg score if both score1 and score2 has it.
                double delta, totalScore;
                if (double.IsNaN(score1) || double.IsNaN(score2))
                {
                    if(useScGs) return; //if a scoreResult has been set using SCGs, dont bother
                    score2 = double.IsNaN(score2) ? GetIcoScoreDynamic(cluster2) : score2;
                    var commonCo = CommonCoCalculation(cluster1, cluster2, _coMatrix);
                    delta = commonCo - Math.Max(score1Ico, score2) ;
                    if (delta < 0) return;
                    totalScore = commonCo;
                }
                else
                {
                    //calculate score first, as is typically faster and might fail more often
                    var comboScore = _evaluator.ScoreCluster(cluster1.Union(cluster2));
                    delta = comboScore - Math.Max(score1, score2);
                    if (delta <= 0) return;

                    var commonCo = CommonCoCalculation(cluster1, cluster2, _coMatrix);
                    totalScore = (1+commonCo) * (comboScore);
                    localScgUsage = true;
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
                            useScGs = useScGs || localScgUsage; //update whether or not to compare using onlu SCGs, used for speedup
                        }
                    }
                }
            });
        
        return scoreResult;
    }


    private List<string[]> SplitAllNegativeClusters(IEnumerable<string[]> clusters)
    {
        Console.WriteLine("Splitting negative bins...");
        return clusters
            .AsParallel()
            .SelectMany(x => 
                GetClusterScoreDynamic(x) < 0.0 ? SplitBin(x) : new[] { x })
            .ToList();
    }

    private static IEnumerable<string[]> SplitBin(IEnumerable<string> cluster)
    {
        return cluster.Select(x => new []{ x });
    }

    
    private static double CommonCoSum(IReadOnlyList<string> cluster1, IReadOnlyList<string> cluster2,
        IReadOnlyDictionary<CoTuple, double> coMatrix)
    {
        var sum = 0.0d;
        for (var i = 0; i < cluster1.Count; i++)
        {
            for (var j = 0; j < cluster2.Count; j++)
            {
                sum += coMatrix.TryGetValue(new CoTuple(cluster1[i], cluster2[j]), out var value ) ? value : 0.0d;
            }
        }
        return sum;
    }
    private static double CommonCoCalculation(IReadOnlyList<string> cluster1, IReadOnlyList<string> cluster2,
        IReadOnlyDictionary<CoTuple, double> coMatrix)
    {
        if (cluster1.Count == 0 || cluster2.Count == 0) return 0.0d;
        return CommonCoSum(cluster1, cluster2, coMatrix) / (cluster1.Count * cluster2.Count);
    }
    
    public static double InternalCommonCoCalculation(IReadOnlyList<string> cluster1, IReadOnlyList<string> cluster2,
        double score1, double score2, IReadOnlyDictionary<CoTuple, double> coMatrix)
    {
        var overlapSize = cluster1.Count * cluster2.Count;
        if (cluster1.Count == 0 || cluster2.Count == 0) return 0.0d;
        var co = CommonCoCalculation(cluster1, cluster2, coMatrix);
        
        return ( (score1 * cluster1.Count) + (co * overlapSize )  + (score2 * cluster2.Count) )
            / (cluster1.Count + overlapSize + cluster2.Count);
    }
    
}

