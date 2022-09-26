using System.Collections;
using ShellProgressBar;

namespace BinChillingTools;


public static class ProgressBarExtensions
{
    public static string FormatMessage(this IProgressBar progressBar, DateTime startTime = default)
    {
        return $"    Remainder: {FormatTimeRemaining(progressBar.CalculateTimeRemainder(startTime))}";
    }
    
    public static void FormattedTick(this IProgressBar progressBar, DateTime startTime = default)
    {
        progressBar.Tick(progressBar.FormatMessage(startTime));
    }

    public static void Finalize(this IProgressBar progressBar)
    {
        progressBar.Tick(newTickCount: progressBar.MaxTicks, DefaultFinishMsg );
    }

    public static string FormatTimeRemaining(TimeSpan span)
    {
        return $"{span.Hours:00}:{span.Minutes:00}:{span.Seconds:00}";
    }
    
    public static readonly string DefaultStartMsg = $"    Remainder: {FormatTimeRemaining(TimeSpan.MaxValue)}";
    public static readonly string DefaultFinishMsg = $"    Remainder: {FormatTimeRemaining(TimeSpan.Zero)}";
    
    public static ProgressBar Create(int count) => 
        new(count, DefaultStartMsg, new ProgressBarOptions()
        {
            BackgroundCharacter = '-',
            BackgroundColor = ConsoleColor.Black,
            CollapseWhenFinished = false,
            DenseProgressBar = false,
            DisplayTimeInRealTime = false,
            DisableBottomPercentage = false,
            EnableTaskBarProgress = false,
            ForegroundColor = ConsoleColor.Green,
            ForegroundColorDone = ConsoleColor.White,
            ForegroundColorError = ConsoleColor.Red
        });
    
    public static ProgressBarOptions DefaultChildOptions() => new()
    {
        BackgroundCharacter = '-',
        BackgroundColor = ConsoleColor.Black,
        CollapseWhenFinished = true,
        DenseProgressBar = false,
        DisplayTimeInRealTime = false,
        DisableBottomPercentage = false,
        EnableTaskBarProgress = false,
        ForegroundColor = ConsoleColor.Green,
        ForegroundColorDone = ConsoleColor.White,
        ForegroundColorError = ConsoleColor.Red
    };


    public static ProgressBar CreateProgressBar<T>(this IReadOnlyCollection<T> collection) => Create(collection.Count);

    public static TimeSpan CalculateTimeRemainder(this IProgressBar progressBar, DateTime startTime)
    {
        var tick = progressBar.CurrentTick > 0 ? progressBar.CurrentTick : 1;
        startTime = startTime < DateTime.Now ? startTime : DateTime.Now;
        var timespan = DateTime.Now.Subtract(startTime).Ticks;
        return TimeSpan.FromTicks((timespan / tick) * (progressBar.MaxTicks - tick));
    }
}

public sealed class ProgressBarEnumerable<T> : IEnumerable<T>
{
    private readonly IEnumerable<T> _enumerable;
    private readonly int _size;
    private readonly Func<int, IProgressBar> _progressBarGenerator;
    private IProgressBar? _bar;

    public ProgressBarEnumerable(IEnumerable<T> enumerable, int size, 
        Func<int, IProgressBar> progressBarGenerator, IProgressBar? bar = null)
    {
        _enumerable = enumerable;
        _size = size;
        _progressBarGenerator = progressBarGenerator;
        _bar = bar;
    }

    private IProgressBar GetProgressBar()
    {
        var bar = _bar;
        if (bar is null) return _progressBarGenerator.Invoke(_size);
        _bar = null;
        return bar;
    }
    
    
    public IEnumerator<T> GetEnumerator()
    {
        return _enumerable is ProgressBarEnumerable<T> 
            ? _enumerable.GetEnumerator() 
            : new ProgressBarEnumerator<T>(_enumerable, GetProgressBar());
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
    
    
}

public sealed class ProgressBarEnumerator<T> : IEnumerator<T>
{
    private readonly IEnumerable<T> _enumerable;
    private readonly IProgressBar _progressBar;
    
    private IEnumerator<T>? _enumerator;
    private DateTime? _startTime;

    public ProgressBarEnumerator(IEnumerable<T> enumerable, IProgressBar progressBar)
    {
        _enumerable = enumerable;
        _progressBar = progressBar;
    }

    public bool MoveNext()
    {
        if (_startTime.HasValue is false)
        {
            _startTime = DateTime.Now;
        }
        
        if (_enumerator is null)
        {
            _enumerator = _enumerable.GetEnumerator();
            
        }
        
        var next = _enumerator.MoveNext();
        if (next) _progressBar.FormattedTick(_startTime.Value);
        else _progressBar.Finalize();

        return next;
    }

    public void Reset()
    {
        _enumerator?.Dispose();
        _enumerator = null;
        _startTime = null;
        _progressBar.Tick(newTickCount:0, _progressBar.FormatMessage(DateTime.Now));
    }

    public T Current => _enumerator is null ? default! : _enumerator.Current!;

    object IEnumerator.Current => Current ?? new object();

    public void Dispose()
    {
        _enumerator?.Dispose();
        _progressBar.Dispose();
    }
}