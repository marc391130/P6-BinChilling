using System.Text;

namespace BinChillingTools;


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


    public IReadOnlyDictionary<CoTuple, double> ReadCoMatrix()
    {
        Console.WriteLine("Reading co matrix...");
        return File.ReadLines(_coPath)
            .Select(x => x.TrimEnd('\n').Split('\t'))
            .Where(x => x.Length == 3)
            .ToDictionary(
                k => new CoTuple(k[0], k[1]),
                v => double.TryParse(v[2], out var val)
                    ? val
                    : throw new FileLoadException("I cannot read value: " + v[2]));
    }


    public IReadOnlyList<string[]> ReadPartition()
    {
        Console.WriteLine("Reading Partitions...");
        return ReadPartitionInternal(_partitionPath).Select(x => x.ToArray()).ToArray();
    }

    private static IEnumerable<ISet<string>> ReadPartitionInternal(string path)
    {
        var group = new Dictionary<string, ISet<string>>();
        
        foreach (var line in File.ReadLines(path)
                     .Select(x => x.TrimEnd('\n').Split('\t'))
                     .Where(x => x.Length == 2))
        {
            if (group.ContainsKey(line[0]))
                group[line[0]].Add(line[1]);
            else
                group[line[0]] = new HashSet<string>() { line[1] };
        }

        return group.Values;
    }

    public static IReadOnlyDictionary<CoTuple, double> ReadSourcePartitionAndBuildCoMatrix(
        IReadOnlyList<string[]> partition, string folder)
    {
        Console.WriteLine("Reading source partitions to build co-matrix");
        var dictionary = new Dictionary<CoTuple, int>();
        var sourcePartitions = new List<IReadOnlyList<ISet<string>>>();
        var totalBinsCount = partition.Sum(x => x.Length);


        foreach (var source in Directory.GetFiles(folder))
        {
            sourcePartitions.Add(ReadPartitionInternal(source).ToList());
        }
        
        using (var progressBar = new ConsoleProgressBar((uint)totalBinsCount))
        {
            foreach (var item1 in partition.SelectMany(static x => x).UseProgressBar(progressBar))
            {
                foreach (var sourcePartition in sourcePartitions)
                {
                    var cluster = sourcePartition
                                      .FirstOrDefault(x => x.Contains(item1))
                                  ?? new HashSet<string>(0);
                    foreach (var item2 in cluster)
                    {
                        var tup = new CoTuple(item1, item2);
                        dictionary[tup] = 1 + (dictionary.TryGetValue(tup, out var v) ? v : 0);
                    }
                }   
            }
        }
        
        var partitionCount = Directory.GetFiles(folder).Length;
        return dictionary.ToDictionary(
            k => k.Key,
            v => ((double)v.Value) / ((double)partitionCount)
        );
    }

    //key is contig
    //value is list of SCGs
    public IReadOnlyDictionary<string, string[]> ReadScgs()
    {
        Console.WriteLine("Reading Scgs...");
        var group = new Dictionary<string, List<string>>();
        foreach (var line in File.ReadLines(_scgPath)
                     .Select(x => x.TrimEnd('\n').Split('\t'))
                     .Where(x => x.Length == 2))
        {
            if (group.ContainsKey(line[0]))
                group[line[0]].Add(line[1]);
            else
                group.Add(line[0], new List<string>(){ line[1] });
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