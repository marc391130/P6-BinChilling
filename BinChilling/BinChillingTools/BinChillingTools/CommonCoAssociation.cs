using System.Collections.ObjectModel;

namespace BinChillingTools;

public static class CommonCoAssociation
{

    private static ICollection<ICollection<int>>? _sharedHashes;
    
    public static IEnumerable<RecalcResult> Build(ICollection<ICollection<int>> clusterLst, ICollection<ICollection<int>> newClusters)
    {
        try
        {
            _sharedHashes = clusterLst;

            var oldNewResults = newClusters.AsParallel()
                .Select(static (x, i) => NewOldCommonCoCalculation(i, x))
                .ToArray();

            _sharedHashes = newClusters;
            
            var newNewResults = newClusters.AsParallel()
                .Select(static (x, i) => NewNewCommonCoCalculation(i, x))
                .ToArray();

            return oldNewResults.Concat(newNewResults);
        }
        finally
        {
            _sharedHashes = null;
        }
    }

    private static RecalcResult NewOldCommonCoCalculation(int index, ICollection<int> cluster)
    {
        if (_sharedHashes is null) throw new ArgumentNullException(nameof(_sharedHashes));

        var result = new SubResult[_sharedHashes.Count];

        for (var i = 0; i < _sharedHashes.Count; i++)
        {
            var clusterHash = _sharedHashes.ElementAt(i);
            if(clusterHash.Count == 0) continue;
            var entry = CommonCoCalculation(cluster, clusterHash);

            result[i] = new SubResult(i, entry);
        }

        return new RecalcResult(index, result);
    }
    
    private static RecalcResult NewNewCommonCoCalculation(int index, ICollection<int> cluster)
    {
        if (_sharedHashes is null) throw new ArgumentNullException(nameof(_sharedHashes));

        var result = new SubResult[_sharedHashes.Count];

        for (var i = index; i < _sharedHashes.Count; i++)
        {
            var clusterHash = _sharedHashes.ElementAt(i);
            if(clusterHash.Count == 0) continue;
            var entry = CommonCoCalculation(cluster, clusterHash);

            result[i] = new SubResult(i, entry);
        }

        return new RecalcResult(index, result);
    }

    private static double CommonCoCalculation(ICollection<int> cluster1, ICollection<int> cluster2)
    {
        var value = (
            from i in cluster1 
            from j in cluster2 select CoMatrix.CombineHash(i, j) 
            into index 
            select CoMatrix.AsReadonly.TryGetValue(index, out var val) ? val : 0.0d)
            .Sum();

        return value / (cluster1.Count + cluster2.Count);
    }
}