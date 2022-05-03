# Auther: Narucs, Dak, Spand, & Trudsliv
# Created: 23-02-2022
# Intended purpose: To get data from fasta files

from io import TextIOWrapper
from random import randrange, random, seed

from AdaptiveEnsembler import AdaptiveClusterEnsembler, Ensembler
from AdaptiveEnsemblerExtensions import MergeRegulator, target_bin_3_4th_count_estimator, QualityMeasuerer
from AdaptiveEnsemblerDomainExtensions import  MergeSCGEvaluator
from BinEvaluator import BinEvaluator
from Binner2 import Binner2
from Cluster import Cluster, Partition, PartitionSet
from typing import Callable
from tqdm import tqdm
from EvaluationMethods import ARIEvaluator, NMIEvaluator
from PartitionSetReader import PartitionSetReader
from ContigReader import ContigFilter, ContigReader
from Domain import ContigData
import Constants as Const 
import argparse
import sys
import os
import scipy as sp
from BinChilling import BinChillingEnsembler, Binner, Chiller, MyLogger

# import numpy
# arr = numpy.zeros((3,3))
# app = numpy.ones(shape=(3,1))
# print("> arr\n", arr)
# print("> app\n", app)
# print(numpy.append(arr , app, axis=1 ))
# raise Exception()

def print_result(file_path: str, parititon: Partition[ContigData]):
    with open(file_path, 'w') as file:
        cluster_lst = list(parititon.values())
        for cluster_idx in range(len(cluster_lst)):
            for item in cluster_lst[cluster_idx]:
                file.write(f"{cluster_idx+1}\t{item.name}\n")


def run(logger: MyLogger, a1:float, a1_min: float, target_cluster_est: int or Callable[[PartitionSet], int],\
        fasta_filepath: str, depth_filepath: str, scg_filepath: str,\
    numpy_cachepath: str, partition_folder: str, output_path: str,\
        max_processors: int or None, chunksize: int, min_contig_len: int, use_old: bool):
    
    #Load data
    contigFilter = ContigFilter(min_contig_len)
    contigReader = ContigReader(fasta_filepath, depth_filepath, scg_filepath, SCG_db_path='../Dataset/Bacteria.ms',\
        numpy_file=numpy_cachepath, max_threads=max_processors)
    partitionSetReader = PartitionSetReader(partition_folder, contigReader, lambda x: x.endswith(".tsv"),\
        contig_filter=contigFilter)
    partition_set = partitionSetReader.read_file()
    
    
    #set up ensembler
    bin_evaluator = BinEvaluator(contigReader.read_total_SCGs_set())
    print('> XD ', len(bin_evaluator.all_SCGs))
    regulator = MergeRegulator(a1_min)  if use_old else\
                MergeSCGEvaluator(a1_min, bin_evaluator, debug=True)
    chiller = Chiller(a1_min, a1, regulator, 0.02, logger)
    binner = Binner2(bin_evaluator, QualityMeasuerer(), 0.75, logger)
    ensembler = BinChillingEnsembler(chiller, binner, bin_evaluator, target_cluster_est, chunksize, max_processors, logger)

    output = ensembler.ensemble(partition_set)
    # #TODO MOVE THIS TO A BETTER PLACE LATER!
   
    # # regulator = MergeRegulator(ensembler.aplha1_min, ensembler.taget_clusters_est)
    # ensembler.merge_regulator = regulator

    # output = ensembler.ensemble(partition_set)
        
    #Display output
    print_result(output_path, output)
    print("Completed successfully")
    all_elements_len = len(partition_set.get_all_elements())
    true_cluster = partition_set[0]
    print(ARIEvaluator.evaluate(output, true_cluster, all_elements_len))
    print(NMIEvaluator.evaluate(output, true_cluster, all_elements_len))
    print(f'Can now be run on Evaluator: it has {all_elements_len} objects!')

    sys.exit(0)
    
def run_old(a1_min, fasta_filepath: str, depth_filepath: str, scg_filepath: str, numpy_cachepath: str, partition_folder: str,\
    output_path: str, threads, chunksize, logfile: TextIOWrapper or None, min_contig_len: int = 0):

    ensembler = AdaptiveClusterEnsembler(
        initial_alpha1_thredshold=1,\
        initial_delta_aplha=0.02,\
        alpha1_min=a1_min,\
        alpha2=0.8,\
        taget_clusters_est=target_bin_3_4th_count_estimator,\
        logfile=logfile,\
        should_log=True,\
        threads=threads,\
        chunksize=chunksize)
    
    contigFilter = ContigFilter(min_contig_len)
    contigReader = ContigReader(fasta_filepath, depth_filepath, scg_filepath, SCG_db_path='../Dataset/Bacteria.ms',\
        numpy_file=numpy_cachepath, max_threads=threads)
    partitionSetReader = PartitionSetReader(partition_folder, contigReader, lambda x: x.endswith(".tsv"),\
        contig_filter=contigFilter)
    partition_set = partitionSetReader.read_file()
    
    regulator = MergeRegulator(ensembler.aplha1_min)
    ensembler.merge_regulator = regulator

    output = ensembler.ensemble(partition_set)
        
    print_result(output_path, output)
    print("Completed successfully")
    all_elements_len = len(partition_set.get_all_elements())
    true_cluster = partition_set[0]
    print(ARIEvaluator.evaluate(output, true_cluster, all_elements_len))
    print(NMIEvaluator.evaluate(output, true_cluster, all_elements_len))
    print(f'Can now be run on Evaluator: it has {all_elements_len} objects!')

    sys.exit(0)

def main():    
    seed(2) # most random number used for seed, chosen by committee
    
    parser = argparse.ArgumentParser(
        prog='ACE',
        description="""ACE BINNING ENSEMBLER""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s WRITE THIS LATER PLZ",
        add_help=True
        )
    
    
    p_args = parser.add_argument_group(title='Contig input (required)', description=None)
    p_args.add_argument('--fasta', metavar='', required=True,\
        dest='fasta', help='path to fasta file of contigs')
    p_args.add_argument('--SCG', metavar='', required=True, \
        dest='SCG', help='Path to single copy genes file (required)')
    p_args.add_argument('--jgi', metavar='', required=False, default=None, \
        dest='JGI', help='path to depth file (either this or --NPZ required)')
    p_args.add_argument('--NPZ', metavar='', required=False, default=None, \
        dest='NPZ', help='path to abundance file in npz format (either this or --jgi required)')
    
    IO_args = parser.add_argument_group(title='IO (required)', description=None)
    IO_args.add_argument('--cache', metavar='', type=str, required=False,\
        dest='numpy_cache', help='Path for cache. If no cache exists at the path, a cache file will be created (highly encuraged)')
    IO_args.add_argument('-P', metavar='', type=str, required=True,\
        dest='partition_folder', help='Path to folder of partition files (Each partition should be its own .tsv file)')
    IO_args.add_argument('--outfile', metavar='', type=str, required=True,\
        dest='outdir', help='Output file of the result (Overrides file if already exists)')
    IO_args.add_argument('-l', metavar='', type=str, required=False, default=None,\
        dest='logdest', help='file to output logfile [default = no log file]')
    IO_args.add_argument('-t', metavar='', type=int, required=False, default=None,\
        dest='threads', help='Max number of threads to use [default <= 8]')
    
    
    ensemble_args = parser.add_argument_group(title='Ensemble variables', description=None)
    ensemble_args.add_argument('-a1', type=float, dest='a1', metavar='',
        default=1, help='initial a1 value, for merging similar clusters [default = 1]')

    ensemble_args.add_argument('-a1min', type=float, dest='a1_min', metavar='',
        default=0.85, help='the minimum threshold for merging clusters. (recommended between 0.5 and 0.9) [default = 0.85]')
    ensemble_args.add_argument('-a2', type=float, dest='a2', metavar='',\
        default=0.85 , help='initial a2 value, for labeling item certainty [default = 0.85]')
    ensemble_args.add_argument('-m', type=int, dest='min_contigs', metavar='',
        default=100, help='minimum contig length to use [default = 100]')
    
    ensemble_args.add_argument('-k', type=int, dest='target_clusters', metavar='',
        default=None, help='The number of bins to target during the process [default = 3rd quartile average]')
    ensemble_args.add_argument('-c', type=int, dest='chunksize', metavar='',
        default=50, help='The chinksize to split a list into when multithreading [default = 50, ignored if -t = 1]')
    ensemble_args.add_argument('--old', type=bool, dest='use_old', metavar='', default=False, \
        help='Use default ACE criterea for breaking merging process')
    
    if(len(sys.argv) <= 1):
        parser.print_help()
        sys.exit()
        
    args = parser.parse_args()
    print(args.use_old)
    ###### BINNING ARGS ######
    
    #fasta file exist 
    fasta_path = os.path.abspath(args.fasta)    
    if os.path.isfile(fasta_path) is False:
        raise FileNotFoundError(fasta_path)
    
    #jgi exist
    abundance_path = None
    if args.JGI is not None:
        abundance_path = os.path.abspath(args.JGI)
        if os.path.isfile(abundance_path) is False:
            raise FileNotFoundError(abundance_path)
    elif args.NPZ is not None:
        abundance_path = os.path.abspath(args.NPZ)
    else:
        raise argparse.ArgumentError(None, 'Either JGI or NPZ option required')
    
    if os.path.isfile(abundance_path) is False:
        raise FileNotFoundError(abundance_path)
    
    #Single copy genes file
    SCG_path = os.path.abspath(args.SCG)
    if os.path.isfile(SCG_path) is False:
        raise FileNotFoundError(SCG_path)
    
    ###### IO ARGS ######
    numpy_cache = None
    if args.numpy_cache is not None:
        numpy_cache: str = os.path.abspath(args.numpy_cache)
        if not numpy_cache.endswith('.npy'):
            raise argparse.ArgumentError(args.numpy_cache, message='Provided cache file does not end in .npy')
    
    #partition folder
    partition_folder = args.partition_folder if args.partition_folder.endswith('/') else args.partition_folder + '/'
    partition_folder = os.path.dirname(partition_folder)
    if partition_folder and os.path.isdir(partition_folder) is False:
        raise NotADirectoryError(partition_folder)
    
    #output file
    outfile: str = os.path.abspath(args.outdir)
    outfile = outfile if outfile.endswith('.tsv') else outfile + '.tsv'
    if os.path.isfile(outfile):
        print(f"Output file '{outfile}' already exists, overriding file when process completes...")
    
    ###### ENSEMBLER ARGS ######
    
    if  0 > args.a1 or args.a1 > 2: 
        raise argparse.ArgumentError(args.a1, "a1 is not in range 0 to 1")
    
    if 0 > args.a1_min or args.a1_min > 1:
        raise argparse.ArgumentError(args.a1_min, "a1_min is not in range 0 to 1") 
    
    if 0 > args.a2 or args.a2 > 1:
        raise argparse.ArgumentError(args.a2, "a2 is not in range 0 to 1") 

    target_clusters = args.target_clusters if args.target_clusters is not None\
        else target_bin_3_4th_count_estimator
    
    if args.chunksize < 1:
        raise argparse.ArgumentError(args.chunksize, 'chunksize must be larger than 0')
    
    try: 
        logfile = open(args.logdest, 'w') if args.logdest is not None else None
        logger = MyLogger(console_log=True, logfile=logfile)

        
                
        # ensembler = AdaptiveClusterEnsembler(
        #         initial_alpha1_thredshold=args.a1,
        #         initial_delta_aplha=0.02,
        #         alpha1_min=args.a1_min,
        #         alpha2=args.a2,
        #         taget_clusters_est=target_clusters,
        #         logfile=logfile,
        #         should_log=True,
        #         threads=args.threads,
        #         chunksize=args.chunksize
        #     )
        if args.use_old:
            run_old(args.a1_min, fasta_path, abundance_path, SCG_path, numpy_cache, partition_folder, 
                    outfile, args.threads, args.chunksize, logfile, args.min_contigs)
        else:
            run(logger, args.a1, args.a1_min, target_clusters, fasta_path, abundance_path, SCG_path,\
                numpy_cache, partition_folder, outfile, args.threads, args.chunksize, args.min_contigs, args.use_old)
    finally:
        if logfile is not None:
            logfile.close()
    
if __name__ == '__main__':
    main()
    

    # # ensembler = AdaptiveClusterEnsembler( \
    # #     initial_alpha1_thredshold=0.6, \
    # #     initial_delta_aplha=0.1, \
    # #     alpha1_min=0.9, \
    # #     alpha2=0.6)
    


    # use_real_data = True
    # candidate_clusters = None
    # if use_real_data:
    #     contigReader = ContigReader(Const.FASTA_FILEPATH, Const.ABUNDANCE_FILEPATH, Const.SCG_FILEPATH, "../Dataset/contigs_numpy.npy")
    #     partitionSetReader = PartitionSetReader("../Dataset/ClusterData/", contigReader, lambda x: x.endswith(".tsv"))
    #     partition_set = partitionSetReader.read_file()

    #     candidate_clusters = ensembler.ensemble(partition_set)
    # else:
    #     def dummy_data(partitions = 2, clusters_per_partition = 100, elements_in_data=1000) -> PartitionSet:
    #         data = [Contig(str(random())) for x in tqdm(range(elements_in_data))]

    #         partition_set = PartitionSet(data)
    #         for x in range(partitions):
    #             partition_set.create_partition()

    #         for partition in range(len(partition_set)):
    #             cluster_max = max(int(randrange(int(clusters_per_partition * 0.1), clusters_per_partition)), 2)
                
    #             for x in data:
    #                 cluster = int(randrange(0,clusters_per_partition))
    #                 partition_set[partition].add(str(cluster), x)
    #         return partition_set

    #     candidate_clusters = ensembler.ensemble(dummy_data(partitions=10, clusters_per_partition=2, elements_in_data=10))

    # print_result("../Dataset/ACE_output.tsv", candidate_clusters)
    # # ensembler.print_to_file("../Dataset/ACE_output.tsv", candidate_clusters)


