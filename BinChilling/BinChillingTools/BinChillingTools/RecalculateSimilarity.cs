using System.Collections.Concurrent;
using System.Security.Cryptography;

namespace BinChillingTools;

public record RecalcResult(int Id, SubResult[] Values);
public record SubResult(int Id2, double Value); 


public static class RecalculateSimilarity
{

    private static IReadOnlyCollection<IReadOnlyCollection<int>>? _clusterHashes;

    
    public static RecalcResult[] Recalculate(
        IEnumerable<int> itemLst, 
        int[][] clusterHashes, 
        Func<int, int[]> clusterListDelegate)
    {
        try
        {
            _clusterHashes = clusterHashes;

            return itemLst
                .AsParallel()
                .Select((hash, index) => 
                    PartialRecalculateSimilarity(index, hash, clusterListDelegate.Invoke(hash)))
                .ToArray();
        }
        finally
        {
            _clusterHashes = null;
        }
    }
    
    private static RecalcResult PartialRecalculateSimilarity(int id, int itemHash, ICollection<int> clusterIndexes)
    {
        var result = new SubResult[clusterIndexes.Count];
        if (_clusterHashes is null) throw new ArgumentNullException(nameof(_clusterHashes));

        for (var i = 0; i < clusterIndexes.Count; i++)
        {
            var cIndex = clusterIndexes.ElementAt(i);
            var value = 0.0;
            var clusterHashes = _clusterHashes.ElementAt(cIndex);
            if(clusterHashes.Count == 0) continue;

            foreach (var item2Hash in clusterHashes)
            {
                var index = CoMatrix.CombineHash(itemHash, item2Hash);
                value += CoMatrix.AsReadonly.TryGetValue(index, out var co) ? co : 0.0d;
            }

            value /= clusterHashes.Count;
            result[i] = new SubResult(cIndex, value);
        }

        return new RecalcResult(id, result);
    }


    
    
    
}
