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

        if (bag.IsEmpty) return 1.0d;
        return bag.Average();
    }
    
    public static IEnumerable<KeyValuePair<T, IReadOnlyList<T>>> Segment<T>(this IReadOnlyList<T> list)
    {
        for (var i = 0; i < list.Count; i++)
        {
            var offset = i + 1;
            var remainder = list.Count - offset;
            yield return new KeyValuePair<T, IReadOnlyList<T>>(
                list[i], 
                new ListSegment<T>(list, offset, remainder)
                );
        }
    }
    
    
    public static IEnumerable<KeyValuePair<T, IReadOnlyList<T>>> ExtendSegment<T>(
        this IEnumerable<KeyValuePair<T, IReadOnlyList<T>>> segment,
        IReadOnlyList<T> extension)
    {
        return segment.Select(p => new KeyValuePair<T, IReadOnlyList<T>>(
            p.Key, new ConcatCollection<T>(p.Value, extension)));
    }
}