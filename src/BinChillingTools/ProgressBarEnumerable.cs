using System.Collections;

namespace BinChillingTools;

public sealed class ProgressBarEnumerable<T> : IEnumerable<T>
{
    private readonly IEnumerable<T> _enumerable;
    private readonly ConsoleProgressBar _progressBar;

    public ProgressBarEnumerable(IEnumerable<T> enumerable, ConsoleProgressBar progressBar)
    {
        _enumerable = enumerable;
        _progressBar = progressBar;
    }

    public IEnumerator<T> GetEnumerator()
    {
        return new ProgressBarEnumerator<T>(_enumerable.GetEnumerator(), _progressBar);
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
}

public sealed class ProgressBarEnumerator<T> : IEnumerator<T>
{
    private readonly ConsoleProgressBar _progressBar;
    private readonly IEnumerator<T> _enumerator;

    public ProgressBarEnumerator(IEnumerator<T> enumerator, ConsoleProgressBar progressBar)
    {
        _enumerator = enumerator;
        _progressBar = progressBar;
    }
    
    public bool MoveNext()
    {
        var next = _enumerator.MoveNext();
        _progressBar.Update();
        return next;
    }

    public void Reset()
    {
        _enumerator.Reset();
        _progressBar.Reset();
    }

    public T Current => _enumerator.Current;

    object IEnumerator.Current => Current!;

    public void Dispose()
    {
        _enumerator.Dispose();
    }
} 