# Auther: Narucs, Dak, Spand, & Trudsliv
# Created: 23-02-2022
# Intended purpose: To get data from fasta files

from io import TextIOWrapper
from random import seed
from typing import List, Dict, Tuple
from functools import partialmethod
from multiprocessing import cpu_count

from AdaptiveEnsembler import AdaptiveClusterEnsembler
from BinChillingBinner import run_binner
from EnsemblerTools import BinLogger, MergeRegulator, target_bin_3_4th_count_estimator, print_result
from BinEvaluator import BinEvaluator
from BinRefiner import *
from ClusterDomain import PartitionSet
from typing import Callable
from tqdm import tqdm
from BinReaders import ContigFilter, SCGReader, PartitionSetReader, ContigReader
import argparse
import sys
import os
from BinChillingEnsembler import BinChillingEnsembler, Chiller, Binner
from time import time
import numpy as np

# import numpy
# arr = numpy.zeros((3,3))
# app = numpy.ones(shape=(3,1))
# print("> arr\n", arr)
# print("> app\n", app)
# print(numpy.append(arr , app, axis=1 ))
# raise Exception()

__global_disable_tqdm = False
tqdm.__init__ = partialmethod(tqdm.__init__, disable=__global_disable_tqdm)


def run_ensemble(logger: BinLogger, a1:float, a1_min: float, target_cluster_est: int or Callable[[PartitionSet], int],\
        fasta_filepath: str, depth_filepath: str, scg_filepath: List[str], db_path: List[str],\
    numpy_cachepath: str, partition_folder: str, output_path: str,\
        max_processors: int or None, chunksize: int, min_contig_len: int):
    
    start_time = time()
    
    #Load data
    scg_reader = SCGReader(scg_filepath, db_path, logger=logger)
    contigFilter = ContigFilter(min_contig_len)
    contigReader = ContigReader(fasta_filepath, scg_reader,\
        depth_file=depth_filepath, numpy_file=numpy_cachepath, max_threads=max_processors, logger=logger)
    partitionSetReader = PartitionSetReader(partition_folder, contigReader, lambda x: x.endswith(".tsv"),\
        contig_filter=contigFilter, logger=logger)
    gamma = partitionSetReader.read_file()
    
    
    #set up ensembler
    bin_evaluator = BinEvaluator(scg_reader.read_MS_scgs())
    # print('>SCG coutn ', len(bin_evaluator.all_SCGs))
    regulator = MergeRegulator(a1_min) if True else\
                MergeSCGEvaluator(a1_min, bin_evaluator, debug=True)
    chiller = Chiller(a1_min, a1, regulator, 0.02, logger)
    binner = Binner(bin_evaluator, chunksize, 0.75, logger)
    ensembler = BinChillingEnsembler(chiller, binner, bin_evaluator, target_cluster_est, chunksize, max_processors, logger)

    output = ensembler.ensemble(gamma)

    print_result(output_path, output)
    
    logger.log("Starting refinement process")
    
    # bin_refiner = BinRefiner(bin_evaluator, 1 / len(gamma), chunksize, logger)
    external_binrefiner = ExternalBinRefiner(len(gamma), output_path + '.ref.tsv', output_path, path.join(os.getcwd(), "comatrix.tmp") , logger)
    external_binrefiner.refine(output)
    
    #Display output
    logger.log(f"Finished bin ensemblement in time {(time() - start_time):0.02f}s")
    
    logger.log("Completed successfully")
    sys.exit(0)
    
def run_ACE(a1_min, fasta_filepath: str, depth_filepath: str, scg_filepath: List[str], db_file: List[str], numpy_cachepath: str, partition_folder: str,\
    output_path: str, threads, chunksize, logfile: TextIOWrapper or None, min_contig_len: int = 0):

    scg_reader = SCGReader(scg_filepath, db_file)   
    contigFilter = ContigFilter(min_contig_len)
    contigReader = ContigReader(fasta_filepath, scg_reader, depth_filepath, scg_filepath,\
        numpy_file=numpy_cachepath, max_threads=threads)
    partitionSetReader = PartitionSetReader(partition_folder, contigReader, lambda x: x.endswith(".tsv"),\
        contig_filter=contigFilter)
    gamma = partitionSetReader.read_file()
    
    ensembler = AdaptiveClusterEnsembler(
        initial_alpha1_thredshold=1,\
        initial_delta_aplha=0.02,\
        alpha1_min=a1_min,\
        alpha2=(1 - (1/len(gamma))),\
        taget_clusters_est=target_bin_3_4th_count_estimator,\
        logfile=logfile,\
        should_log=True,\
        threads=threads,\
        chunksize=chunksize)
    
    
    regulator = MergeRegulator(ensembler.aplha1_min)
    ensembler.merge_regulator = regulator

    output = ensembler.ensemble(gamma)
        
    print_result(output_path, output)
    print("Completed successfully")
    sys.exit(0)

def main():    
    parser = argparse.ArgumentParser(
        prog='BinChilling',
        description="""BinChilling Ensemble and Binning algorithm implementation. 
        Includes an implementation of the generic ACE ensembler (althought it only works for binning, for now)""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s WRITE THIS LATER PLZ",
        add_help=True
        )
    
    req_args = parser.add_argument_group(title='Required args', description=None)
    req_args.add_argument('--Module', '-M', metavar='', required=False, type=str, default='Ensemble',\
        dest='module', choices=['Bin', 'Ensemble', 'ACE'], help="The module to run ('Ensemble', 'Bin', 'ACE'). (REQUIRED)"
        )
    
    p_args = parser.add_argument_group(title='Common IO options (used by all modules)', description=None)
    p_args.add_argument('--fasta' ,'-f', metavar='', type=str, required=True,\
        dest='fasta', help='path to fasta file of contigs [required]')
    p_args.add_argument('--SCG', '-g', nargs='+', metavar='', required=True, \
        dest='SCG', help='Path to single copy marker genes files (.tsv format) [requires at least 1]')
    p_args.add_argument('--MS', '-MS', dest='MS', nargs='+', default=None, \
        help='Path to MS files of genes. If no ms files supplied, the union of all SCGs are used [optional]', metavar='', required=False)
    p_args.add_argument('--cache', '-c', metavar='', type=str, required=False,\
        dest='numpy_cache', help='Path for contig cache. If no cache exists at the path, a cache file will be created (optional, but highly encuraged)')
    
    p_args.add_argument('--outfile', '-o', metavar='', type=str, required=True,\
        dest='outdir', help='Output file of the result (Overrides file if already exists) [required]')
    p_args.add_argument('--log','-l', metavar='', type=str, required=False, default=None,\
        dest='logdest', help='file to output logfile [default = no log file]')
    
    common_args = parser.add_argument_group(title='Common args (used by all modules)', description=None)
    common_args.add_argument('--threads', '-t', metavar='', type=int, required=False, default=None,\
        dest='threads', help='Max number of threads to use [default = all]')
    common_args.add_argument('-a1', type=float, dest='a1', metavar='',
        default=1.0, help='initial a1 value, for merging similar clusters (1.0 recommended) [default = 1.0]')
    common_args.add_argument('--a1min', '-am', type=float, dest='a1_min', metavar='',
        default=0.90, help='the minimum threshold for merging clusters. (recommended between 0.5 and 0.99) [default = 0.90]')
    common_args.add_argument('--chunksize', '-ch', type=int, dest='chunksize', metavar='',
        default=500, help='The chunksize to split a list into when multithreading [default = 500, ignored if -t = 1]')
    common_args.add_argument('--minlen', '-m', type=int, dest='minlen', metavar='',
        default=1000, help='minimum contig length to use [default = 1000]')
    common_args.add_argument('--rand', '-r', type=int, dest='rand', default=2, \
        help='The seed to use for random number generation. Set seed to 0 to use random seed (default = 2)')
    
    ensemble_args = parser.add_argument_group(title='Ensemble args (only used by ACE and Ensemble modules)', description='Used by both BinChillingEnsembler and ACE')
    ensemble_args.add_argument('--Partitions', '-p', metavar='', type=str, default=None, required=False,\
        dest='partition_folder', help='Path to folder of partition files (Each partition should be its own .tsv file) [required]')
    
    ace_args = parser.add_argument_group(title='ACE args', description='Used by both BinChillingEnsembler and ACE')
    ace_args.add_argument('-a2', type=float, dest='a2', metavar='',\
        default=0.85 , help='initial a2 value, for labeling item certainty [default = 0.85]')
    ace_args.add_argument('-k', type=int, dest='target_clusters', metavar='',
        default=None, help='The number of bins to target during the process [default = 3rd quartile partitions avg]')
    
    binner_args = parser.add_argument_group(title='Binner args')
    binner_args.add_argument('--jgi', '-j', metavar='', required=False, default=None, \
        dest='JGI', help='path to depth file (either this or --NPZ/-z required)')
    binner_args.add_argument('--NPZ', '-z', metavar='', required=False, default=None, \
        dest='NPZ', help='path to abundance file in npz format (either this or --jgi/-j required)')
    binner_args.add_argument('--out_part', '-op', metavar='', required=False, default=None, type=str,\
        help='Folder to output the generated partition before ensemblement [optional]', dest='partition_outdir')
    binner_args.add_argument('--min', '-p0', type=int, dest='min_p', default=3, \
        help='Minimum number of partitions to generate (default = 3)', metavar='', required=False)
    binner_args.add_argument('--max', '-p1', type=int, dest='max_p', default=100, \
        help='Maximum number of partitions to generate (default = 100)', metavar='', required=False)
    binner_args.add_argument('--stepsize', '-s', type=int, dest='stepsize', default=5,\
         help='The amount of bins to add to the next generated partition (default = 5)', metavar='', required=False)
    binner_args.add_argument('--algorithm', '-A',  dest='algorithm', default='PartialSeed', choices=['Kmeans', 'PartialSeed', 'Hierarchical', 'Random'], type=str, metavar='',  required=False, \
        help="Algorithm to use for generating partitions. Choices: 'Kmeans', 'PartialSeed', 'PartialSeed', 'Hierarchical', 'Random'. (default = PartialSeed)\n" + \
             "KMeans: sklearn's implementation using KMeans++ to initialize centers.\n" +\
             "PartialSeed: the KMeans initialization method proposed in the MetaBinner paper.\n" +\
             "Hierarchical: sklearns AgglomerativeClustering module with 'cannot-link' constraints.\n" +\
             "Random: embedding is randomly generated between each partition generation" )
    
    if(len(sys.argv) <= 1):
        parser.print_help()
        sys.exit()
        
    args = parser.parse_args()
    ###### MODULE ARGS ######
    run_binner_mod = args.module == 'Bin'
    run_ACE_mod = args.module == 'ACE'
    run_ensemble_mod = args.module == 'Ensemble'
    if (run_binner_mod or run_ACE_mod or run_ensemble_mod) is False:
        raise ValueError('No module selected. See help for details')
    
    ###### Common IO ARGS ######
    #fasta file exist 
    fasta_path = os.path.abspath(args.fasta)    
    if os.path.isfile(fasta_path) is False:
        raise FileNotFoundError(fasta_path)
    
    #SCG files
    SCG_path = []
    scg = [args.SCG] if isinstance(args.SCG, str) else args.SCG
    if len(scg) == 0:
        raise ValueError('No SCG file supplied')
    for file in scg:
        if os.path.isfile(file) is False:
            raise FileNotFoundError(file)
        SCG_path.append(os.path.abspath(file))
    
    #MS file
    MS_path = []
    if args.MS is not None:
        MS_path = [args.MS] if isinstance(args.MS, str) else args.MS
        for i in range(len(MS_path)):
            MS_path[i] = os.path.abspath(MS_path[i])
            if os.path.isfile(MS_path[i]) is False:
                raise FileNotFoundError(MS_path[i])
    
    #jgi exist
    abundance_path = None
    if run_binner_mod and args.JGI is None and args.NPZ is None:
        raise ValueError('No abundance supplied. Either JGI or NPZ option required for binning')
    elif run_binner_mod:
        if args.JGI is not None:
            abundance_path = os.path.abspath(args.JGI)
            if os.path.isfile(abundance_path) is False:
                raise FileNotFoundError(abundance_path)
        elif args.NPZ is not None:
            abundance_path = os.path.abspath(args.NPZ)

        if os.path.isfile(abundance_path) is False:
            raise FileNotFoundError(abundance_path)
    #endif
    
    numpy_cache = None
    if args.numpy_cache is not None:
        numpy_cache: str = os.path.abspath(args.numpy_cache)
        if not numpy_cache.endswith('.npy'):
            raise ValueError(args.numpy_cache, message='Provided cache file does not end in .npy')
    
    
    #partition folder
    partition_folder = None
    if (run_ensemble_mod or run_ACE_mod) and args.partition_folder is None :
        raise ValueError('folder containing partitions required for bin ensemblement')
    elif (run_ensemble_mod or run_ACE_mod):
        partition_folder = args.partition_folder if args.partition_folder.endswith('/') else args.partition_folder + '/'
        partition_folder = os.path.dirname(partition_folder)
        if partition_folder and os.path.isdir(partition_folder) is False:
            raise NotADirectoryError(partition_folder)
    
    #output file
    outfile: str = os.path.abspath(args.outdir)
    outfile = outfile if outfile.endswith('.tsv') else outfile + '.tsv'
    if os.path.isfile(outfile):
        print(f"Output file '{outfile}' already exists, overriding file when process completes...")
    
    #logfile
    logdest = None
    if args.logdest is not None:
        logdest = os.path.abspath(args.logdest)
        if os.access(os.path.dirname(logdest), os.W_OK) is False:
            raise PermissionError(f'Python doesnt have permission to write logfile to {logdest}')
    
    ###### COMMON ARGS ######
    
    
    threads = min(args.threads, cpu_count()) if args.threads is not None else cpu_count()
    if threads < 1:
        raise ValueError('Cannot run on less than 1 thread.')
        
    
    a1 = args.a1
    if  (0.0 < a1 <= 1.0) is False: 
        raise ValueError(f"a1 of {a1} is not in range 0.0 to 1.0")
    
    a1min = args.a1_min
    if (0 <= a1min <= 1.0) is False:
        raise ValueError(f"a1_min of {a1min} is not in range 0 to 1") 
    
    a2 = args.a2
    if (0 <= a2 <= 1) is False:
        raise ValueError(f"a2 of {a2} is not in range 0 to 1") 


    target_clusters = args.target_clusters if args.target_clusters is not None\
        else target_bin_3_4th_count_estimator
    if callable(target_clusters) is False and target_clusters < 1:
        raise ValueError('Cannot have target clusters be less than 1')
    
    chunksize = args.chunksize
    if chunksize < 1:
        raise ValueError('chunksize must be larger than 0')
    
    minlen = int(args.minlen)
    if minlen < 0:
        raise ValueError('minlen cannot be negative')
    
    if args.rand != 0:
        np.random.seed(args.rand)
        seed(args.rand)
    
    
    ###### Binner ARGS ######
    partition_outdir = args.partition_outdir if args.partition_outdir is not None else None 
    if partition_outdir and os.path.isdir(partition_outdir) is False:
            raise NotADirectoryError(partition_outdir)
    
    if args.min_p < 1:
        raise ValueError(f'minimum partition generation (min / p0) should be at least 1. recieved {args.min_p}')
    if args.min_p > args.max_p: 
        raise ValueError("min nr. of partitions is higher than max")
    
    if args.stepsize < 1:
        raise ValueError(f'binner stepsize shoud be at least 1. recieved {args.stepsize}')
    
    if args.algorithm not in ['Kmeans', 'PartialSeed', 'Hierarchical', 'Random']:
        raise ValueError('Binner partition generator is not recognized. See help for details')
    
    try: 
        logfile = open(logdest, 'w') if logdest is not None else None
        logger = BinLogger(console_log=True, logfile=logfile)


        if run_binner_mod:
            run_binner(a1min, args.min_p, args.max_p, minlen, args.stepsize, args.algorithm, fasta_path, abundance_path,\
                SCG_path, MS_path, outfile, numpy_cache, chunksize,  logger, 300, partition_outdir)
        elif run_ensemble_mod:
            run_ensemble(logger, args.a1, args.a1_min, target_clusters, fasta_path, abundance_path, scg, MS_path,\
                numpy_cache, partition_folder, outfile, args.threads, args.chunksize, minlen)
        elif run_ACE_mod:
            run_ACE(args.a1_min, fasta_path, abundance_path, SCG_path, MS_path, numpy_cache, partition_folder, 
                    outfile, args.threads, args.chunksize, logfile, minlen)
        else:
            raise Exception('module not recognized')
            
    finally:
        if logfile is not None:
            logfile.close()
    
if __name__ == '__main__':
    main()