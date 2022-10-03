using System.Collections;
using ShellProgressBar;

namespace BinChillingTools;

public class ProgressBarEnumerable<T> : IEnumerable<T>
{
    private readonly IEnumerable<T> _enumerable;
    private readonly IProgressBar _progressBar;

    public ProgressBarEnumerable(IEnumerable<T> enumerable, IProgressBar progressBar)
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

public class ProgressBarEnumerator<T> : IEnumerator<T>
{
    private readonly IProgressBar _progressBar;
    private readonly IEnumerator<T> _enumerator;
    private DateTime _startTime;

    public ProgressBarEnumerator(IEnumerator<T> enumerator, IProgressBar progressBar)
    {
        _enumerator = enumerator;
        _progressBar = progressBar;
    }
    
    public bool MoveNext()
    {
        var next = _enumerator.MoveNext();
        _progressBar.FormattedTick(_startTime);
        return next;
    }

    public void Reset()
    {
        _enumerator.Reset();
        _startTime = DateTime.Now;
        _progressBar.Tick(newTickCount:0, ProgressBarHandler.DefaultStartMsg);
    }

    public T Current => _enumerator.Current;

    object IEnumerator.Current => Current!;

    public void Dispose()
    {
        _enumerator.Dispose();
    }
} 