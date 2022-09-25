using ShellProgressBar;

namespace BinChillingTools;

public static class EnumerableCoExtensions
{
    public static IEnumerable<T> UseProgressBar<T>(this IEnumerable<T> enumerable, 
        int size, string message = "", ProgressBarOptions? options = null)
    {
        return new ProgressBarEnumerable<T>(enumerable, size, message, options);
    }

    public static IEnumerable<KeyValuePair<T, IEnumerable<T>>> Segment<T>(this T[] list, Comparison<T>? comparison = null)
    {
        if (comparison is not null) Array.Sort(list, comparison);

        for (var i = 0; i < list.Length; i++)
        {
            var offset = i + 1;
            var remainder = list.Length - offset;
            yield return new KeyValuePair<T, IEnumerable<T>>(
                list[i], 
                new ArraySegment<T>(list, offset, remainder)
                );
        }
    }
    
    public static IEnumerable<KeyValuePair<T, IEnumerable<T>>> SmartGroupSegment<T>(this IReadOnlyList<T> list, IReadOnlyList<T> segment,
        Comparison<T>? comparison = null)
    {
        for (var i = 0; i < list.Count; i++)
        {
            var selfSegment = list.Skip(i+1);
            yield return new KeyValuePair<T, IEnumerable<T>>(list[i], selfSegment.Concat(segment));
        }
    }
}