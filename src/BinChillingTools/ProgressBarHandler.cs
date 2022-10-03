using System.Collections;
using ShellProgressBar;

namespace BinChillingTools;


public static class ProgressBarHandler
{
    public static IEnumerable<T> UseProgressBar<T>(this IEnumerable<T> enumerable, IProgressBar progressBar)
    {
        return new ProgressBarEnumerable<T>(enumerable, progressBar);
    }

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
    
    public static readonly string DefaultStartMsg =  $"    Remainder: {FormatTimeRemaining(TimeSpan.MaxValue)}";
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
    
    public static readonly ProgressBarOptions DefaultChildOptions = new()
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
