using System.Runtime.InteropServices;
using RGiesecke.DllExport;

namespace BinChillingTools;

public static class PythonInterface
{

    [DllExport(nameof(BuildCoMatrix), CallingConvention.Cdecl)]
    public static void BuildCoMatrix(IEnumerable<CoEntry> input, dynamic source)
    {
        CoMatrix.Build(input);
    }

    [DllExport(nameof(ClearCoMatrix), CallingConvention.Cdecl)]
    public static void ClearCoMatrix()
    {
        CoMatrix.Clear();
    }

    [DllExport(nameof(RecalculateSim), CallingConvention.Cdecl)]
    public static IEnumerable<RecalcResult> RecalculateSim(IEnumerable<int> itemLst, 
        int[][] clusterHashes, 
        Func<int, int[]> clusterListDelegate)
    {
        return RecalculateSimilarity.Recalculate(itemLst, clusterHashes, clusterListDelegate);
    }
    
    [DllExport(nameof(RecalculateCommonCo), CallingConvention.Cdecl)]
    public static IEnumerable<RecalcResult> RecalculateCommonCo(ICollection<ICollection<int>> clusterLst, ICollection<ICollection<int>> newClusters)
    {
        return CommonCoAssociation.Build(clusterLst, newClusters);
    }
    
}