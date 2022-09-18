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
        Console.WriteLine($"totalScgs: {_totalScgs}");
    }


    public double ScoreCluster(IEnumerable<string> cluster)
    {
        var map   = CalculateScgMap(cluster);
        var completeness    = 100.0d * Completeness(map);
        var contamination   = 100.0d * Contamination(map);
        var megabinPenality = 100.0d * MegabinPenality(map);

        var score = completeness - Math.Pow(contamination, 2) - Math.Sqrt(megabinPenality);
        return score;
    }


    private Dictionary<string, int> CalculateScgMap(IEnumerable<string> cluster)
    {
        var dic = new Dictionary<string, int>();
        foreach (var scg in cluster.SelectMany(x => _scgDict.TryGetValue(x, out var lst) 
                     ? lst 
                     : Array.Empty<string>()))
        {
            if (dic.TryGetValue(scg, out var value))
            {
                dic[scg] = value + 1;
            }
            else
            {
                dic.Add(scg, 1);
                // dic[scg] = 1;
            }
        }
        return dic;
    }

    private double Completeness(Dictionary<string, int> scgMap)
    {
        float counter = scgMap.Count; //nr. of unique scgs in cluster
        float divisor = _totalScgs; //nr. of scgs total

        return divisor != 0 ? counter / divisor : 0;
    }
    
    private double Contamination(Dictionary<string, int> scgMap)
    {
        float counter = scgMap.Count(x => x.Value > 1);
        float divisor = scgMap.Count;

        return divisor != 0 ? counter / divisor : 0;
    }
    
    private double MegabinPenality(Dictionary<string, int> scgMap)
    {
        float counter = scgMap.Values.Sum() - scgMap.Count;
        float divisor = _totalScgs;

        return divisor != 0 ? counter / divisor : 0;
    }
    
}