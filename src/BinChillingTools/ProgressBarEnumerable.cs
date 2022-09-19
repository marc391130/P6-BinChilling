using System.Collections;
using ShellProgressBar;

namespace BinChillingTools;

public static class EnumerableProgressBarExtensions
{
    public static IEnumerable<T> UseProgressBar<T>(this IEnumerable<T> enumerable, 
        int size, string message = "", ProgressBarOptions? options = null)
    {
        return new ProgressBarEnumerable<T>(enumerable, size, message, options);
    }
}

public sealed class ProgressBarEnumerable<T> : IEnumerable<T>
{
    private readonly IEnumerable<T> _enumerable;
    private readonly int _size;
    private readonly string _message;
    private readonly ProgressBarOptions? _options;

    public ProgressBarEnumerable(IEnumerable<T> enumerable, int size, string message = "", ProgressBarOptions? options = null)
    {
        _enumerable = enumerable;
        _size = size;
        _message = message;
        _options = options;
    }
    public IEnumerator<T> GetEnumerator()
    {
        return _enumerable is ProgressBarEnumerable<T> 
            ? _enumerable.GetEnumerator() 
            : new ProgressBarEnumerator<T>(_enumerable, _size, _message, _options);
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
    
    
}

public sealed class ProgressBarEnumerator<T> : IEnumerator<T>
{
    private const ConsoleColor DefaultColor = ConsoleColor.White;
    private readonly IEnumerable<T> _enumerable;
    private readonly int _size;
    private readonly string _message;
    private readonly ProgressBarOptions? _options;
    private ProgressBar _bar;
    private IEnumerator<T>? _enumerator;

    public ProgressBarEnumerator(IEnumerable<T> enumerable, int size, string message, ProgressBarOptions? options = null)
    {
        _enumerable = enumerable;
        _size = size;
        _message = message;
        _options = options;
        _bar = options is null 
            ? new ProgressBar(size, message, DefaultColor) 
            : new ProgressBar(size, message, options);
    }

    private static ProgressBar BuildProgressBar(int size, string message, ProgressBarOptions? options = null) => options is null 
        ? new ProgressBar(size, message, DefaultColor) 
        : new ProgressBar(size, message, options);
        
    public bool MoveNext()
    {
        if (_enumerator is null)
        {
            _enumerator = _enumerable.GetEnumerator();
        }
        
        _bar.Tick();
        return _enumerator.MoveNext();
    }

    public void Reset()
    {
        _enumerator = null;
        _bar.Dispose();
        _bar = BuildProgressBar(_size, _message, _options);
    }

    public T Current => _enumerator is null ? default : _enumerator.Current!;

    object IEnumerator.Current => Current ?? new object();

    public void Dispose()
    {
        _enumerator?.Dispose();
        _bar.Dispose();
    }
}