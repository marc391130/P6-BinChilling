namespace BinChillingTools;

public readonly struct CoTuple
{
    private readonly string _item1;
    private readonly string _item2;
    
    
    public CoTuple(string item1, string item2)
    {
        if (string.Compare(item1, item2, StringComparison.InvariantCulture) <= 0)
        {
            _item1 = item1;
            _item2 = item2;
        }
        else
        {
            _item2 = item1;
            _item1 = item2;
        }
    }

    public override int GetHashCode()
    {
        return HashCode.Combine(_item1, _item2);
    }

    public override bool Equals(object? obj)
    {
        if (obj is not CoTuple tup) return false;

        return _item1 == tup._item1 && _item2 == tup._item2;
    }

    public static bool operator ==(CoTuple left, CoTuple right)
    {
        return left.Equals(right);
    }

    public static bool operator !=(CoTuple left, CoTuple right)
    {
        return !(left == right);
    }
}