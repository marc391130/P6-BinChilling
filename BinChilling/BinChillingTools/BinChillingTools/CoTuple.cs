namespace BinChillingTools;

public readonly struct CoTuple<T> where T : notnull
{
    public readonly T Item1;
    public readonly T Item2;

    public CoTuple(T item1, T item2)
    {
        var i1 = item1.GetHashCode();
        var i2 = item2.GetHashCode();

        if (i1 < i2)
        {
            Item1 = item1;
            Item2 = item2;
        }
        else
        {
            Item1 = item2;
            Item2 = item1;
        }
    }

    public override int GetHashCode()
    {
        var i1 = Item1.GetHashCode();
        var i2 = Item2.GetHashCode();

        return HashCode.Combine(i1, i2);
    }

    public override bool Equals(object? obj)
    {
        if (obj is CoTuple<T> tup) 
            return tup.Item1.Equals(Item1) &&
                   tup.Item2.Equals(Item2);

        return false;
    }

    public static bool operator ==(CoTuple<T> left, CoTuple<T> right)
    {
        return left.Equals(right);
    }

    public static bool operator !=(CoTuple<T> left, CoTuple<T> right)
    {
        return !(left == right);
    }
}