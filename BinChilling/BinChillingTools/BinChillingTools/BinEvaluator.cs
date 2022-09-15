namespace BinChillingTools;

public sealed class BinEvaluator
{
    private readonly IReadOnlyDictionary<string, string[]> _scgDict;
    private readonly int _totalScgs;

    public BinEvaluator(IReadOnlyDictionary<string, string[]> scgDict)
    {
        _scgDict = scgDict;
        _totalScgs = scgDict.Values
            .SelectMany(static x => x)
            .Distinct()
            .Count();
    }


    public double ScoreCluster(IEnumerable<string> cluster)
    {
        var map   = CalculateScgMap(cluster);
        var completeness    = 100 * Completeness(map);
        var contamination   = 100 * Contamination(map);
        var megabinPenality = 100 * MegabinPenality(map);

        return completeness - Math.Pow(contamination, 2) - Math.Sqrt(megabinPenality);
    }


    private Dictionary<string, int> CalculateScgMap(IEnumerable<string> cluster)
    {
        var dic = new Dictionary<string, int>();
        foreach (var scg in cluster.SelectMany(x => _scgDict.TryGetValue(x, out var lst) ? lst : Array.Empty<string>()))
        {
            if (dic.TryGetValue(scg, out var value))
            {
                dic[scg] = value + 1;
            }
            else
            {
                dic[scg] = 1;
            }
        }
        return dic;
    }

    private double Completeness(Dictionary<string, int> scgMap)
    {
        var counter = scgMap.Count; //nr. of unique scgs in cluster
        var divisor = _totalScgs; //nr. of scgs total

        return divisor != 0 ? counter / divisor : 0;
    }
    
    private double Contamination(Dictionary<string, int> scgMap)
    {
        var counter = scgMap.Count(x => x.Value > 1);
        var divisor = scgMap.Count;

        return divisor != 0 ? counter / divisor : 0;
    }
    
    private double MegabinPenality(Dictionary<string, int> scgMap)
    {
        var counter = scgMap.Values.Sum() - scgMap.Count;
        var divisor = _totalScgs;

        return divisor != 0 ? counter / divisor : 0;
    }
    
}