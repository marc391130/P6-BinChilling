using ShellProgressBar;

namespace BinChillingTools;

public static class EnumerableCoExtensions
{
    public static IEnumerable<T> UseProgressBar<T>(this IReadOnlyCollection<T> collection)
    {
        return new ProgressBarEnumerable<T>(collection, collection.Count, ProgressBarExtensions.Create);
    }

    /// <summary>
    /// ASSUMES progressBar.MaxTicks is the size of the collection
    /// </summary>
    /// <param name="enumerable"></param>
    /// <param name="progressBar"></param>
    /// <typeparam name="T"></typeparam>
    /// <returns></returns>
    public static IEnumerable<T> UseProgressBar<T>(this IEnumerable<T> enumerable, IProgressBar progressBar)
    {
        return new ProgressBarEnumerable<T>(enumerable, progressBar.MaxTicks, ProgressBarExtensions.Create, progressBar);
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
    
    public static IEnumerable<KeyValuePair<T, IReadOnlyCollection<T>>> SmartGroupSegment<T>(this IReadOnlyList<T> list,
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