using System.Collections.Concurrent;

namespace BinChillingTools;

public static class EnumerableCoExtensions
{
    public static double InternalCoScore(this IReadOnlyList<string[]> clusters,
        IReadOnlyDictionary<CoTuple, double> coMatrix, ConsoleProgressBar? progressBar = null)
    {
        var sum = 0.0d;
        if (clusters.Count <= 0) return sum;
        for (int k = 0; k < clusters.Count; k++)
        {
            progressBar?.Update();
            var cluster = clusters[k];
            if (cluster.Length == 0) continue;
            var localSum = 0.0d;
            var localCount = 0;
            for (int i = 0; i < cluster.Length; i++)
            {
                for (int j = i+1; j < cluster.Length; j++)
                {
                    localCount++;
                    localSum += coMatrix.TryGetValue(new CoTuple(cluster[i], cluster[j]), out var value) ? value : 0.0d;
                }
            }
            if(localCount <= 0) continue;
            sum += localSum / localCount;
        }
        return sum / clusters.Count;
    }

    public static double InternalCo(this IReadOnlyList<string> cluster, IReadOnlyDictionary<CoTuple, double> coMatrix)
    {
        var sum = 0.0d;
        var count = 0;
        for (int i = 0; i < cluster.Count; i++)
        {
            for (int j = i+1; j < cluster.Count; j++)
            {
                count++;
                sum += coMatrix.TryGetValue(new CoTuple(cluster[i], cluster[j]), out var value) ? value : 0.0d;
            }
        }

        if (count <= 0) return 0.0d;
        return sum / ((double)count);
    }
    
    public static double InternalCoScoreParallel(this IEnumerable<string[]> clusters,
        IReadOnlyDictionary<CoTuple, double> coMatrix, ConsoleProgressBar? progressBar = null)
    {
        var bag = new ConcurrentBag<double>();
        var enumerable = progressBar is null 
            ? clusters 
            : clusters.UseProgressBar(progressBar);

        Parallel.ForEach(enumerable.Where(x => x.Length > 0),
            (cluster) => bag.Add(InternalCo(cluster, coMatrix)));

        return bag.Average();
    }
    
    public static IEnumerable<KeyValuePair<T, IReadOnlyCollection<T>>> Segment<T>(this T[] list, Comparison<T>? comparison = null)
    {
        if (comparison is not null) Array.Sort(list, comparison);

        for (var i = 0; i < list.Length; i++)
        {
            var offset = i + 1;
            var remainder = list.Length - offset;
            yield return new KeyValuePair<T, IReadOnlyCollection<T>>(
                list[i], 
                new ArraySegment<T>(list, offset, remainder)
                );
        }
    }
    
    public static IEnumerable<KeyValuePair<T, IReadOnlyCollection<T>>> SmartGroupSegment<T>(this List<T> list,
        IReadOnlyList<T> segment, Comparison<T>? comparison = null)
    {
        for (var i = 0; i < list.Count; i++)
        {
            var offset = i + 1;
            var remainder = list.Count - offset;
            var selfSegment = new ListSegment<T>(list, offset, remainder);
            
            yield return new KeyValuePair<T, IReadOnlyCollection<T>>(
                list[i], 
                new ConcatCollection<T>(selfSegment, segment)
                );
        }
    }
}