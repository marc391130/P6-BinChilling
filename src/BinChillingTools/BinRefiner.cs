using System.Collections.Concurrent;
using ShellProgressBar;

namespace BinChillingTools;


public readonly record struct CoResult(double Score, string[] Cluster);

public sealed class BinRefiner
{
    private readonly double _minCoValue;
    private readonly IReadOnlyDictionary<CoTuple, double> _coMatrix;
    private readonly BinEvaluator _evaluator;
    private readonly ConcurrentDictionary<int, double> _scoreDict = new();

    public BinRefiner(double minCoValue, IReadOnlyDictionary<CoTuple, double> coMatrix, BinEvaluator evaluator)
    {
        _minCoValue = minCoValue;
        _coMatrix = coMatrix;
        _evaluator = evaluator;
    }


    public IReadOnlyList<IReadOnlyList<string>> Refine(IReadOnlyDictionary<string, string[]> partition)
    {
        Console.WriteLine($"Starting refinement using {partition.Count} clusters");
        
        _scoreDict.Clear();
        var namedClusters = partition
            .Select(x => x.Value)
            .ToList();

        var clusterList = namedClusters as IReadOnlyList<string[]>;
        var removeSet = new HashSet<int>();


        while (true)
        {
            var newClusters = RefineClustersMultiThreaded_NoBreak(clusterList, removeSet);
            Console.WriteLine($"Refined {removeSet.Count} clusters into {newClusters.Count} new clusters.");
            if (newClusters.Count < 2) 
                break;
            

            foreach (var clusterHash in removeSet)
            {
                _scoreDict.Remove(clusterHash, out _);
            }
            clusterList = clusterList
                .Where(x => removeSet.Contains(x.GetHashCode()) is false)
                .Concat(newClusters)
                .ToList();
            
            removeSet.Clear();
            GC.Collect();
        }
        
        _scoreDict.Clear();
        return clusterList;
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


    private IReadOnlyList<string[]> RefineClusters(
        IReadOnlyList<string[]> clusters, ISet<int> skipSet)
    {
        var newClusters = new List<string[]>();
        for (var i = 0; i < clusters.Count; i++)
        {
            var cluster1 = clusters[i];
            if(skipSet.Contains(cluster1.GetHashCode()) || cluster1.Length == 0) continue;
            var score1 = GetClusterScoreDynamic(cluster1);
            
            for (var j = i+1; j < clusters.Count; j++)
            {
                var cluster2 = clusters[j];
                if(skipSet.Contains(cluster2.GetHashCode()) || cluster2.Length == 0) continue;
                
                if(CommonCoCalculation(cluster1, cluster2, _coMatrix) <= _minCoValue) continue;
                var score2 = GetClusterScoreDynamic(cluster2);
                
                var comboScore = _evaluator.ScoreCluster(cluster1.Concat(cluster2));
                if (comboScore > Math.Max(score1, score2))
                {
                    // Console.WriteLine("Combo: " + comboScore + " | c1: " + score1 + " | c2: " + score2);
                    
                    skipSet.Add(cluster1.GetHashCode());
                    skipSet.Add(cluster2.GetHashCode());
                    newClusters.Add( cluster1.Concat(cluster2).ToArray() );
                    break;
                }
            }
            
            if ( skipSet.Contains(cluster1.GetHashCode())) continue;
            if (score1 < 0.0d)
            {
                skipSet.Add(cluster1.GetHashCode());
                newClusters.AddRange( SplitBin( cluster1 ) );
                continue;
            }
            
            // newClusters.Add(cluster1);
        }

        return newClusters;
    }
    
    
    private IReadOnlyList<string[]> RefineClustersMultiThreaded(
        IReadOnlyList<string[]> clusters, ISet<int> skipSet)
    {
        var newClusters = new List<string[]>();
        for (var i = 0; i < clusters.Count; i++)
        {
            var cluster1 = clusters[i];
            if(skipSet.Contains(cluster1.GetHashCode()) || cluster1.Length == 0) continue;
            var score1 = GetClusterScoreDynamic(cluster1);

            Parallel.For(i + 1, clusters.Count, (j, state) =>
            {
                var cluster2 = clusters[j];
                if(skipSet.Contains(cluster2.GetHashCode()) || cluster2.Length == 0) return;

                if(CommonCoCalculation(cluster1, cluster2, _coMatrix) <= _minCoValue) return;
                var score2 = GetClusterScoreDynamic(cluster2);
                
                var comboScore = _evaluator.ScoreCluster(cluster1.Concat(cluster2));
                if (comboScore > Math.Max(score1, score2)) //not locked, as to have maximum throughput
                {
                    lock (skipSet) //lock skipset to go in here.
                    {
                        if(state.IsStopped) return; //assure we are the first to find match
                        
                        skipSet.Add(cluster1.GetHashCode());
                        skipSet.Add(cluster2.GetHashCode());
                        newClusters.Add( cluster1.Concat(cluster2).ToArray() );
                        
                        state.Break(); //break loop
                        return;
                    }
                }
            });

            if ( skipSet.Contains(cluster1.GetHashCode())) continue;
            if (score1 < 0.0d)
            {
                skipSet.Add(cluster1.GetHashCode());
                newClusters.AddRange( SplitBin( cluster1 ) );
                continue;
            }
            
            // newClusters.Add(cluster1);
        }

        return newClusters;
    }
    
    private static ProgressBar CreateProgressBar(int count) => 
        new(count, "Refining clusters...", new ProgressBarOptions()
    {
        BackgroundCharacter = '-',
        BackgroundColor = ConsoleColor.Black,
        CollapseWhenFinished = false,
        DenseProgressBar = true,
        DisplayTimeInRealTime = true,
        DisableBottomPercentage = true,
        EnableTaskBarProgress = false,
        ForegroundColor = ConsoleColor.Green,
        ForegroundColorDone = ConsoleColor.White,
        ForegroundColorError = ConsoleColor.Red
    });

    private static TimeSpan CalculateRemainder(DateTime startTime, ProgressBarBase progressBar)
    {
        startTime = startTime < DateTime.Now ? startTime : DateTime.Now;
        return TimeSpan.FromTicks(DateTime.Now.Subtract(startTime).Ticks / (progressBar.CurrentTick) 
                                  * progressBar.MaxTicks);
    }

    private IReadOnlyList<string[]> RefineClustersMultiThreaded_NoBreak(
        IReadOnlyList<string[]> clusters, ISet<int> skipSet)
    {
        using var progressBar = CreateProgressBar(clusters.Count);
        var startTime = DateTime.Now;
        var newClusters = new List<string[]>();

        for (var i = 0; i < clusters.Count; i++)
        {
            progressBar.Tick($"{progressBar.Percentage:0.00}%|{progressBar.CalculateTimeRemainder(startTime)}");
            var cluster1 = clusters[i];
            if(skipSet.Contains(cluster1.GetHashCode()) || cluster1.Length == 0) continue;
            var score1 = GetClusterScoreDynamic(cluster1);

            CoResult? scoreResult = null;
            Parallel.For(i + 1, clusters.Count, (j) =>
            {
                var cluster2 = clusters[j];
                if(skipSet.Contains(cluster2.GetHashCode()) || cluster2.Length == 0) return;
                
                var commonCo = CommonCoCalculation(cluster1, cluster2, _coMatrix);
                if(commonCo <= _minCoValue) return;
                var score2 = GetClusterScoreDynamic(cluster2);
                
                var comboScore = _evaluator.ScoreCluster(cluster1.Concat(cluster2));
                var delta = comboScore - Math.Max(score1, score2);
                if (delta > 0 && (scoreResult is null || comboScore > scoreResult.Value.Score)) //not locked, as to have maximum throughput
                {
                    lock (skipSet) //lock skipset to go in here.
                    {
                        scoreResult = new CoResult(comboScore*commonCo, cluster2);
                    }
                }
            });
            
            //Add best cluster, if any
            if (scoreResult.HasValue)
            {
                skipSet.Add(cluster1.GetHashCode());
                skipSet.Add(scoreResult.Value.GetHashCode());
                newClusters.Add( cluster1.Concat(scoreResult.Value.Cluster).ToArray() );
            }

            if ( skipSet.Contains(cluster1.GetHashCode())) continue;
            if (score1 < 0.0d)
            {
                skipSet.Add(cluster1.GetHashCode());
                newClusters.AddRange( SplitBin( cluster1 ) );
                continue;
            }
            
            // newClusters.Add(cluster1);
        }
        
        progressBar.Tick(newTickCount: progressBar.MaxTicks, 
            $"{progressBar.Percentage:0.00}% | { CalculateRemainder(startTime, progressBar) }");
        return newClusters;
    }

    private static IEnumerable<string[]> SplitBin(IEnumerable<string> cluster)
    {
        return cluster.Select(x => new []{ x });
    }

    private static double CommonCoCalculation(string[] cluster1, string[] cluster2,
        IReadOnlyDictionary<CoTuple, double> coMatrix)
    {
        var sum = 0.0d;
        if (cluster1.Length == 0 || cluster2.Length == 0) return sum;
        for (int i = 0; i < cluster1.Length; i++)
        {
            for (int j = 0; j < cluster2.Length; j++)
            {
                sum += coMatrix.TryGetValue(new CoTuple(cluster1[i], cluster2[j]), out var value ) ? value : 0.0d;
            }
        }

        return sum / (cluster1.Length * cluster2.Length);
    }
    
}

