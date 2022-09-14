namespace BinChillingTools;

public record CoEntry(int Item1Hash, int Item2Hash, double Value);

public static class CoMatrix
{
    private static readonly Dictionary<int, double> CoMatrixWrite = new Dictionary<int, double>();
    public static IReadOnlyDictionary<int, double> AsReadonly => CoMatrixWrite;

    
    public static void Build(IEnumerable<CoEntry> input)
    {
        foreach (var (item1, item2, value) in input)
        {
            CoMatrixWrite.Add(HashCode.Combine(item1, item2), value);
        }
    }

    public static void Clear()
    {
        CoMatrixWrite.Clear();
    }
    
    public static int CombineHash(int h1, int h2)
    {
        var small = Math.Min(h1, h2);
        var big = Math.Max(h1, h2);

        return HashCode.Combine(small, big);
    }
}