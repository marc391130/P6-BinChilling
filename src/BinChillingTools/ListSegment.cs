using System.Collections;

namespace BinChillingTools;

public readonly struct ListSegment<T> : IReadOnlyList<T>, ICollection<T>
{
    private readonly IReadOnlyList<T> _list;
    private readonly int _offset;
    private readonly int _count;

    public ListSegment(IReadOnlyList<T> list, int offset, int count = -1)
    {
        _list = list ?? throw new NullReferenceException(nameof(list));
        _offset = offset;
        _count = count >= 0 ? count : list.Count - offset;
        if ((uint)offset > (uint)_list.Count || (uint)count > (uint)(list.Count - offset))
            throw new ArgumentOutOfRangeException($"list of size {_list.Count} cannot have offset {_offset} and count {_count}");
    }
    
    public IEnumerator<T> GetEnumerator()
    {
        for (var i = _offset; i < _offset + _count; i++)
        {
            yield return _list[i];
        }
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }

    public void Add(T item) => throw new NotSupportedException();

    public void Clear() => throw new NotSupportedException();

    public bool Contains(T item)
    {
        using var it = GetEnumerator();
        return Enumerable.Contains(this, item);
    }

    public void CopyTo(T[] array, int arrayIndex)
    {
        var i = 0;
        foreach (var item in this)
        {
            array[arrayIndex + i] = item;
            i++;
        }
    }

    public bool Remove(T item) => throw new NotSupportedException();

    public int Count => _count;
    public bool IsReadOnly => true;
    
    private T ElementAt(int index)
    {
        if (_count > 0 && index < _count)
        {
            return _list[_offset + index];
        }

        throw new ArgumentOutOfRangeException(nameof(index));
    }


    public T this[int index] => ElementAt(index);
}