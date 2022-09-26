using System.Collections;

namespace BinChillingTools;

public class ConcatCollection<T> : IReadOnlyCollection<T>
{
    private readonly IReadOnlyCollection<T> _l1;
    private readonly IReadOnlyCollection<T> _l2;

    public ConcatCollection(IReadOnlyCollection<T> l1, IReadOnlyCollection<T> l2)
    {
        _l1 = l1;
        _l2 = l2;
    }
    
    public IEnumerator<T> GetEnumerator()
    {
        return _l1.Concat(_l2).GetEnumerator();
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public bool Contains(T item)
    {
        return _l1.Contains(item) || _l2.Contains(item);
    }

    public int Count => _l1.Count + _l2.Count;
}