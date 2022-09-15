using System.Text;

namespace BinChillingTools;


public static class CoMatrix
{
    public static int CoHash(string item1, string item2)
    {
        var i1 = item1.GetHashCode();
        var i2 = item2.GetHashCode();
        
        return i1 < i2 
            ? HashCode.Combine(i1, i2)
            : HashCode.Combine(i2, i1); 
    }
}

public sealed class BinChillingFileReader
{
    private readonly string _coPath;
    private readonly string _partitionPath;
    private readonly string _scgPath;

    
    
    public BinChillingFileReader(
        string coPath, 
        string partitionPath, 
        string scgPath)
    {
        _coPath = coPath;
        _partitionPath = partitionPath;
        _scgPath = scgPath;
    }


    public IReadOnlyDictionary<int, double> ReadCoMatrix()
    {
        var co = File.ReadLines(_coPath)
            .Select(x => x.TrimEnd('\n').Split('\t'))
            .Where(x => x.Length == 3)
            .ToDictionary(
                k => CoMatrix.CoHash(k[0], k[1]),
                v => double.TryParse(v[2], out var val)
                    ? val
                    : throw new FileLoadException("I cannot fucking read value: " + v[2]));
        
        Console.WriteLine("size: " + co.Count);
        return co;
    }


    public IReadOnlyDictionary<string, string[]> ReadPartition()
    {
        var group = new Dictionary<string, List<string>>();
        foreach (var line in File.ReadLines(_partitionPath)
                     .Select(x => x.TrimEnd('\n').Split('\t'))
                     .Where(x => x.Length == 2))
        {
            if (group.ContainsKey(line[0]))
                group[line[0]].Add(line[1]);
            else
                group[line[0]] = new List<string>() {line[1]};
        }
        return group.ToDictionary(
            k => k.Key, 
            v => v.Value.ToArray());
    }


    //key is contig
    //value is list of SCGs
    public IReadOnlyDictionary<string, string[]> ReadScgs()
    {
        var group = new Dictionary<string, List<string>>();
        foreach (var line in File.ReadLines(_scgPath)
                     .Select(x => x.TrimEnd('\n').Split('\t'))
                     .Where(x => x.Length == 2))
        {
            if (group.ContainsKey(line[0]))
                group[line[0]].Add(line[1]);
            else
                group[line[0]] = new List<string>() {line[1]};
        }
        return group.ToDictionary(
            k => k.Key, 
            v => v.Value.ToArray());
    }


    public void WritePartition(string output, IEnumerable<IEnumerable<string>> clusterLst)
    {
        if(File.Exists(output)) File.Delete(output);

        using var file = File.Create(output);
        using var fileWriter = new StreamWriter(file);

        var i = 1;
        foreach (var cluster in clusterLst)
        {
            foreach (var item in cluster)
            {
                var line = i + "\t" + item;
                fileWriter.WriteLine(line);
            }

            i++;
        }
    }
}