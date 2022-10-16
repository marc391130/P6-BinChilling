﻿namespace BinChillingTools;

public static class Program
{
    private const int NArgs = 6;
    public static int Main(string[] args)
    {
        if (args.Length >= 1 && args[0] == "__test__")
        {
            Console.WriteLine("REFINER MODULE LOADED CORRECTLY");
            Console.WriteLine("__REFINER_MODULE_CONNECTED__"); // do not change, used by python to see if module is callable.
            return 0;
        }
        
        Console.WriteLine($"Args:\n - {string.Join("\n - ", args)}");
        
        if (args.Length < NArgs) throw new ArgumentException(
            "No arguments supplied. Require partition path, co path and scg path, partition counts and output filepath + name, (in this order)");
        if (args.Length > NArgs) Console.WriteLine("More than " + NArgs +" args supplied. Only first " + NArgs +" used");

        var partitionPath = args[0];
        var coPath = args[1];
        var scgPath = args[2];
        // var partitionCount = 1 / double.Parse(args[3]);
        var outputFile = args[4];
        var partitionSizeUpperBound = int.Parse(args[5]);

        var reader = new BinChillingFileReader(coPath, partitionPath, scgPath);

        var evaluator = new BinEvaluator(reader.ReadScgs());
        var coMatrix = reader.ReadCoMatrix();
        
        var refiner = new BinRefiner(partitionSizeUpperBound, coMatrix, evaluator);

        var partition = reader.ReadPartition();
        var refinedPartition = refiner.Refine(partition);
        
        Console.WriteLine($"Finished refinement with {refinedPartition.Count} total clusters.");
        reader.WritePartition(outputFile, refinedPartition);
        
        return 0;
    }

}