

using System.Text;

namespace BinChillingTools;

public class ConsoleProgressBar : IDisposable
{
    private readonly uint _consoleWidth;

    /// <summary>
    /// Total amount of work. Used for calculating current percentage complete .
    /// </summary>
    public uint TotalUnitsOfWork { get; private set; }
    
    public uint CurrentUnitsOfWork { get; private set; }
    
    public DateTime StartTime { get; set; } = DateTime.Now;

    public double Percentage => Math.Clamp(((Math.Min((double)CurrentUnitsOfWork, (double)TotalUnitsOfWork)) 
                                            / ((double)TotalUnitsOfWork)), 0d, 1d );

    public bool Finalized { get; private set; }

    /// <summary>
    /// Constructor.
    /// </summary>
    /// <param name="totalUnitsOfWork">Total amount of work. Used for calculating current percentage complete.</param>
    /// <param name="consoleWidth">the width of the console to draw the progress bar on</param>
    public ConsoleProgressBar(
        uint totalUnitsOfWork,
        uint consoleWidth = 80 )
    {
        TotalUnitsOfWork = totalUnitsOfWork;
        _consoleWidth = consoleWidth;

        Console.WriteLine();
        // _unitsOfWorkPerProgressBlock = (float) TotalUnitsOfWork / _widthInCharacters;
    }

    public void Reset()
    {
        Update(0);
    }
    
    public void Final()
    {
        Update(TotalUnitsOfWork+1);
        Console.WriteLine();
        Finalized = true;
    }
    
    public void Update()
    {
        Update(CurrentUnitsOfWork+1);
    }

    
    private TimeSpan CalculateTimeRemainder()
    {
        return TimeSpan.FromTicks(( (DateTime.Now.Subtract(StartTime)).Ticks / CurrentUnitsOfWork) * (TotalUnitsOfWork - CurrentUnitsOfWork));
    }
    
    /// <summary>
    /// Draws progress bar.
    /// </summary>
    /// <param name="currentUnitOfWork">Current unit of work in relation to TotalUnitsOfWork.</param>
    public void Update(uint currentUnitOfWork)
    {
        CurrentUnitsOfWork = Math.Min(currentUnitOfWork, TotalUnitsOfWork);

        var leftText = $"{100*Percentage:n2}%|";
        var rightText = $"| {CurrentUnitsOfWork}/{TotalUnitsOfWork} " +
                        $"[{CalculateTimeRemainder():hh\\:mm\\:ss}<{DateTime.Now.Subtract(StartTime):hh\\:mm\\:ss}]";
        var textLength = leftText.Length + rightText.Length + 1;    
        
        var completeProgressBlocks = (uint) Math.Round(Percentage * ((float)(_consoleWidth - textLength) ) );

        var completed = string.Empty;
        var nonCompleted = string.Empty;

        if (_consoleWidth - textLength >= 2)
        {
            completed = new string('#', count: (int)completeProgressBlocks);
            nonCompleted = new string('-', count: Math.Max(((int)_consoleWidth - textLength) - ((int)completeProgressBlocks), 0));

        }
        
        // UpdateText(text + completed + nonCompleted);
        Console.Write("\r{0}{1}{2}{3}", leftText, completed, nonCompleted, rightText);
    }
    

    public void Dispose()
    {
        GC.SuppressFinalize(this);
        if (!Finalized) Console.WriteLine();
    }
}