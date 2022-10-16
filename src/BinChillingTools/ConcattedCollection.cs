using System.Collections;

namespace BinChillingTools;

public sealed class ConcatCollection<T> : IReadOnlyList<T>
{
    private readonly IReadOnlyList<T> _l1;
    private readonly IReadOnlyList<T> _l2;

    public ConcatCollection(IReadOnlyList<T> l1, IReadOnlyList<T> l2)
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

    public T this[int index] => index < _l1.Count ? _l1[index] : _l2[index];
}