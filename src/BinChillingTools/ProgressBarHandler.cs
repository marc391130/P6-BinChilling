using System.Collections;

namespace BinChillingTools;


public static class ProgressBarHandler
{
    public static IEnumerable<T> UseProgressBar<T>(this IEnumerable<T> enumerable, ConsoleProgressBar progressBar)
    {
        return new ProgressBarEnumerable<T>(enumerable, progressBar);
    }

    public static ConsoleProgressBar CreateBar(int size)
    {
        return new ConsoleProgressBar((uint)size, (uint)Math.Max(Console.WindowWidth, 1));
    }
}
