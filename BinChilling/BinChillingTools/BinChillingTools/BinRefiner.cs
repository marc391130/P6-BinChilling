using System.Collections.Concurrent;

namespace BinChillingTools;

public sealed class BinRefiner
{
    private readonly double _minCoValue;
    private readonly IReadOnlyDictionary<int, double> _coMatrix;
    private readonly BinEvaluator _evaluator;
    private readonly ConcurrentDictionary<int, double> _scoreDict = new();

    public BinRefiner(double minCoValue, IReadOnlyDictionary<int, double> coMatrix, BinEvaluator evaluator)
    {
        _minCoValue = minCoValue;
        _coMatrix = coMatrix;
        _evaluator = evaluator;
    }


    public IReadOnlyList<IReadOnlyList<string>> Refine(IReadOnlyDictionary<string, string[]> partition)
    {
        _scoreDict.Clear();
        var namedClusters = partition
            .Select(x => x.Value)
            .ToList();

        var newClusters = namedClusters as IReadOnlyList<string[]>;
        var removeSet = new HashSet<int>();


        while (true)
        {
            removeSet.Clear();
            newClusters = RefineClusters(newClusters, removeSet);
            if (removeSet.Count == 0) 
                break;
            
            Console.WriteLine($"Refined {removeSet.Count} clusters into {newClusters.Count} new clusters.");

            foreach (var clusterHash in removeSet)
            {
                _scoreDict.Remove(clusterHash, out _);
            }
        }
        
        _scoreDict.Clear();
        return newClusters;
    }


    private double GetClusterScoreDynamic(IEnumerable<string> cluster)
    {
        var hash = cluster.GetHashCode();
        if(_scoreDict.TryGetValue(hash, out var score)) return score;
        score = _evaluator.ScoreCluster(cluster);
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
                var score2 = GetClusterScoreDynamic(cluster2);
                
                if(CommonCoCalculation(cluster1, cluster2) <= _minCoValue) continue;
                
                var comboScore = _evaluator.ScoreCluster(cluster1.Concat(cluster2));
                if (comboScore >= Math.Max(score1, score2))
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
            
            newClusters.Add(cluster1);
        }

        return newClusters;
    }

    private static IEnumerable<string[]> SplitBin(IEnumerable<string> cluster)
    {
        return cluster.Select(x => new []{ x });
    }

    private double CommonCoCalculation(IReadOnlyCollection<string> cluster1, IReadOnlyCollection<string> cluster2)
    {
        if (cluster1.Count == 0 || cluster2.Count == 0) return 0;
        var sum = 0.0d;
        
        foreach (var item1 in cluster1)
        {
            foreach (var item2 in cluster2)
            {
                sum += _coMatrix.TryGetValue(CoMatrix.CoHash(item1, item2), out var value) ? value : 0.0d;
            }
        }
        
        return sum / (cluster1.Count * cluster2.Count);
    }
    
}

