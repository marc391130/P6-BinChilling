using ShellProgressBar;

namespace BinChillingTools;

public static class EnumerableCoExtensions
{
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