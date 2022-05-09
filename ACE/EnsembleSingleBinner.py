# Auther: Marucs, Bak, Andreas, & Frederik
# Created: 23-02-2022
# Intended purpose: To bin

from io import TextIOWrapper
from sklearn.metrics import silhouette_score
from AdaptiveEnsemblerExtensions import MergeRegulator, QualityMeasuerer, print_result, target_bin_3_4th_count_estimator
from Cluster import Partition, Cluster, PartitionSet
from ContigReader import ContigFilter, ContigReader
import numpy as np
import os
import shutil
import argparse
import random
import Constants as const
from scipy.sparse import lil_matrix, csr_matrix
from sklearn.cluster import KMeans, AgglomerativeClustering
from Domain import ContigData, compute_3_4_scg_count, count_scgs
from typing import List, Dict, TypeVar, Tuple
from tqdm import tqdm
from BinChilling import BinChillingEnsembler, Chiller, MyLogger
from BinEvaluator import BinEvaluator
from Binner2 import Binner2
from BinRefiner import BinRefiner
from time import time
import sys
import Assertions as Assert

CAHCE_DIR = f'./cache_{time()}'

def compute_constraints(data: Dict[int, ContigData]) -> csr_matrix:
    matrix_shape = (len(data), len(data))
    matrix = lil_matrix(np.zeros(shape=matrix_shape), dtype=float, shape=matrix_shape, copy=True)
    lst = list(data.items())
    
    print('Building connectivity matrix...')
    for i1 in tqdm(range(len(lst))):
        index1, contig1 = lst[i1]
        for i2 in range(i1+1, len(lst)):
            index2, contig2 = lst[i2]
            if len(contig1.SCG_genes.intersection(contig2.SCG_genes)) == 0:
                matrix[index1, index2] = 1.0
                
    return matrix.tocsr(copy=True)

def compute_partition_range(features: List[List[float]], weigths: List[float], contigs: List[ContigData], stepsize: int, min_partitions: int or None = None, max_partitions: int or None = None) -> Tuple[int, int]:
    k0 = max(compute_3_4_scg_count(contigs), 2)
    end_k = min(len(features), (k0 + max_partitions if max_partitions is not None else k0*3 + 2*stepsize))
    
    k_val1, k_val2 = np.NINF, np.NINF

    
    data = features
    best_K, best_sil_val = k0+1, np.NINF
    k1, k2 = None, None
    min_partitions = min_partitions+1 if min_partitions is not None else 1
    
    print(f'Searching for best K value. (max value searched will be {end_k})')
    for k in range(k0+1, end_k, stepsize):
        kmeans = KMeans(n_clusters=k, init="k-means++", n_init=10, random_state=7, max_iter=300)
        labels = kmeans.fit_predict(data, sample_weight=weigths)
        score = silhouette_score(data, labels)
        print(f"k value of '{k}'/{end_k} got score: {score}")
        if score > best_sil_val:
            best_sil_val = score
            best_K = k
        elif k > (k0 + min_partitions):
            if k1 is None: 
                k1, k_val1 = k, best_sil_val
                best_sil_val = np.NINF
            elif k2 is None: 
                k2, k_val2 = k, best_sil_val
                break
    
    best_K = k1 if k_val1 > k_val2 else k2
    print(f'Found range partition range {k0} to {best_K}')
    return (k0, best_K)
                                                            
def transform_contigs_to_features(items: List[ContigData]) \
    -> Tuple[Dict[int, ContigData], List[List[float]], List[float], csr_matrix]:
        #returns: contig index map, features, weights, constraint matrix
    index_map, features, weights = {}, [], []
    for i in range(len(items)):
        item = items[i]
        index_map[i] = item
        features.append( item.composition.AsNormalizedFeatureList() )
        weights.append( item.avg_abundance )
    return (index_map, features, weights, compute_constraints(index_map))


def partial_seed_init(scg_count: Dict[str, int], n_clusters: int, index_map: Dict[int, ContigData], featuers: List[List[float]])-> np.ndarray:
    
    scg_sorted, i = list(scg_count.items()), 0
    sorted_contigs: List[Tuple[int, ContigData]] = sorted(index_map.items(), key=lambda x: x[1].contig_length, reverse=True)
    indexes: set[int] = set()
    while True:
        scg, count = scg_sorted[ np.random.randint(0, len(scg_sorted)) ]
        
        for index, contig in sorted_contigs:
            if scg in contig.SCG_genes:
                indexes.add(index)
        if len(indexes) >= n_clusters:
            break
        i += 1
    list_result = [featuers[x] for x in sorted(indexes, key=lambda x: index_map[x].contig_length, reverse=True)[0:n_clusters]] 
    matrix = np.array(list_result, dtype=np.float16)
    return matrix
    

#returns labels
def run_clustering_method(method: str, n_clusters: int, scg_count: Dict[str, int], index_map: Dict[int, ContigData], featuers: List[List[float]], weigths: List[float], constraints: csr_matrix, max_iter: int):
    if method == 'Kmeans':
        return KMeans(n_clusters=n_clusters, max_iter=max_iter, random_state=n_clusters).fit_predict(featuers, sample_weight=weigths)
    if method == 'PartialSeed':
        return KMeans(n_clusters=n_clusters, init=partial_seed_init(scg_count, n_clusters, index_map, featuers), n_init=1, max_iter=max_iter, random_state=0).fit_predict(featuers, sample_weight=weigths)
    if method == 'Hierarchical':
        return AgglomerativeClustering(n_clusters=n_clusters, connectivity=constraints, memory=CAHCE_DIR).fit_predict(featuers)
    if method == 'Random':
        return KMeans(n_clusters=n_clusters, max_iter=max_iter).fit_predict(np.random.rand( len(featuers), 32 ) )
    raise Exception(f'method name {method} not recognized')
    

def run_clustering(partition: Partition, data: Tuple[Dict[int, ContigData], List[List[float]], List[Tuple[int, int]]],\
    method: str, n_clusters: int, scg_count: Dict[str, int] , max_iter: int = 300) -> Tuple[Partition, float]:
    
    index_map, featuers, weights, constraints = data
    labels = run_clustering_method(method, n_clusters, scg_count, index_map, featuers, weights,  constraints, max_iter)
    for i in range(len(labels)):
        partition.add(str(labels[i]), index_map[i])
    return partition

def run(a1min: float, min_partitions_gamma: int, max_partitions_gamma: int, min_contig_len:int, stepsize:int, method: str,\
    fasta_file: str, abundance_file: str, gene_file: List[str], bacteria_file: str, output_file: str, numpy_cache: str, chunksize: int, LList: List[int],\
    logfile: TextIOWrapper or None, partition_outdir: str or None = None):
    
    start_time = time()
    logger = MyLogger(logfile=logfile)
    filter = ContigFilter(min_contig_len)
    contig_reader = ContigReader(fasta_file, abundance_file,\
        gene_file, bacteria_file, numpy_file=numpy_cache, enable_analyse_contig_comp=True)
    
    contigs = [x for x in list(contig_reader.read_file_fast(load_SCGs=True).values()) if filter.predicate(x)] 
    if any( (len(x.SCG_genes) > 0 for x in contigs) ) is False:
        raise Exception('The set of contigs contains no SCGs')
    if any( (x.__has_analysis__() for x in contigs) ) is False:
        raise Exception('The contigs did not contain an analysis of composition. '+\
            'This might be because the cache file is generated using fast analysis. '+\
            'Try deleting cache file, or provide a different cache file path')
    
    gamma = PartitionSet(contigs)
    feature_data = transform_contigs_to_features(contigs)
    
    scg_count = count_scgs(contigs)
    try:
        k_min, k_max = compute_partition_range(feature_data[1], feature_data[2], contigs, stepsize, min_partitions_gamma, max_partitions_gamma)
        k_max = max(k_min + min_partitions_gamma, k_max)
        
        logger.log(f'Generating partitions with {k_min} to {k_max} bins...')
        for k in tqdm(range(k_min, k_max)):
            partition = run_clustering(gamma.create_partition(), feature_data, method, k, scg_count, max_iter=300)
            if partition_outdir is not None:
                print_result(os.path.join(partition_outdir, 'Partition_'+str(k) +'.tsv'), partition)
        #end_loop
    finally:
        shutil.rmtree(os.path.abspath(CAHCE_DIR), ignore_errors=True)
        
    
    bin_evaluator = BinEvaluator(contig_reader.read_total_SCGs_set(), LList)
    bin_refiner = BinRefiner(bin_evaluator, (1.0 / len(gamma)), logger)
    chiller = Chiller(a1min, 1.0, MergeRegulator(a1min), 0.02, logger)
    binner = Binner2(bin_refiner, bin_evaluator, QualityMeasuerer(), logger=logger)
    ensembler = BinChillingEnsembler(chiller, binner, bin_evaluator, chunksize=chunksize, target_clusters_est=target_bin_3_4th_count_estimator, logger=logger)

    final_partition = ensembler.ensemble(gamma)
    print_result(output_file, final_partition)
    
    logger.log(f"Completed binning in time: {(time() - start_time):0.02f}s")

def main():
    parser = argparse.ArgumentParser(
        prog='BinChillingBinner',
        description="""BinChillilng INDIVIDUAL ENSEMBLER""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s WRITE THIS LATER PLZ",
        add_help=True
        )
    
    f_args = parser.add_argument_group('Files/IO')
    f_args.add_argument('--fasta', '-f', type=str, dest='fasta', \
        help='Path to fasta file', metavar='', required=True)
    f_args.add_argument('--SCG', '-S', type=str, dest='SCG', \
        help='Path to SCG file', metavar='', required=True)
    f_args.add_argument('--depth', '-D', type=str, dest='depth', \
        help='Path to Depth file', metavar='', required=True)
    f_args.add_argument('--cache', '-c', metavar='', required=False, dest='numpy_cache', type=str, \
        help='Path for cache (not required). If no cache exists at the path, a cache file will be created (highly encuraged)')
    f_args.add_argument('--output', '-o', type=str, dest='output', \
        help='Path to output file. If file exists, it will be overwritten', metavar='', required=True)
    f_args.add_argument('--MS', '-MS', type=str, dest='MS', \
        help='Path to MS files of genes', metavar='', required=True)
    f_args.add_argument('--out_part', '-op', metavar='', required=False, default=None, type=str,\
        help='Folder to output  partition before ensemblement', dest='partition_outdir')
    f_args.add_argument('--log', '-l', type=str, dest='logfile', default=None, \
        help='Path to output file', metavar='', required=False)
    
    p_args = parser.add_argument_group('Tuning')
    p_args.add_argument('--rand', '-r', type=int, dest='rand', default=2, \
        help='The seed to use for random number generation. Set argument to 0 in order to set seed to specific value (default = 2)')
    p_args.add_argument('--a1min', '-a', type=float, dest='a1min', default=None, \
        help='Minimum a1 value between 0 and 1 to use when merging clusters (default 0.90)', metavar='', required=False)
    p_args.add_argument('--min', '-p0', type=int, dest='min', default=3, \
        help='Minimum number of partitions to use (default = 3)', metavar='', required=False)
    p_args.add_argument('--max', '-p1', type=int, dest='max', default=25, \
        help='Maximum number of partitions to use (default = 25)', metavar='', required=False)
    p_args.add_argument('--chunksize', '-H', type=int, dest='chunksize', default=400, \
        help='Chunksize to use while multiprocessing. Only impacts performance, higher = more memory usage (default=400)', metavar='', required=False)
    p_args.add_argument('--LList', '-L', type=int, nargs='+', dest='LList', default=[1900000, 6500000],\
         help='List of common contig lengths', metavar='', required=False)
    p_args.add_argument('--minSize', '-m', type=int, dest='minSize', default=1000, \
        help='Minimum size of contig to use in binning (default = 1000)', metavar='', required=False),
    p_args.add_argument('--stepsize', '-z', type=int, dest='stepsize', default=5,\
         help='The amount of bins to add to the next generated partition (default = 5)', metavar='', required=False)
    p_args.add_argument('--method', '-A',  dest='method', default=False, choices=['Kmeans', 'PartialSeed', 'Hierarchical', 'Random'], type=str, metavar='',  required=False, \
        help="Algorithm to use for generating partitions. Choices: 'Kmeans', 'PartialSeed', 'PartialSeed', 'Hierarchical', 'Random'. (default = PartialSeed)\n" + \
             "KMeans: sklearn's implementation using KMeans++ to initialize centers.\n" +\
             "PartialSeed: the KMeans initialization method proposed in the MetaBinner paper.\n" +\
             "Hierarchical: sklearns AgglomerativeClustering module with 'cannot-link' constraints.\n" +\
             "Random: embedding is randomly generated between each partition generation" )
    
    if(len(sys.argv) <= 1):
        parser.print_help()
        sys.exit()
    
    args = parser.parse_args()
    
    fasta_path = os.path.abspath(args.fasta)    
    if os.path.isfile(fasta_path) is False:
        raise FileNotFoundError(fasta_path)
    
    #jgi exist
    abundance_path = os.path.abspath(args.depth)  
    if os.path.isfile(abundance_path) is False:
        raise FileNotFoundError(fasta_path)
    
    #Single copy genes file
    SCG_path = os.path.abspath(args.SCG)
    if os.path.isfile(SCG_path) is False:
        raise FileNotFoundError(SCG_path)
    SCG_path = [SCG_path]
    
    MS_path = os.path.abspath(args.MS)
    if os.path.isfile(MS_path) is False:
        raise FileNotFoundError(MS_path)
    MS_path = [MS_path]
    
    outfile: str = os.path.abspath(args.output)
    outfile = outfile if outfile.endswith('.tsv') else outfile + '.tsv'
    if os.path.isfile(outfile):
        print(f"Output file '{outfile}' already exists, overriding file when process completes")
    
    partition_outdir = None
    if args.partition_outdir is not None:
        partition_outdir = args.partition_outdir if args.partition_outdir.endswith('/') else args.partition_outdir + '/'
        partition_outdir = os.path.dirname(partition_outdir)
        if os.path.isdir(partition_outdir) is False:
            raise NotADirectoryError(partition_outdir)
    
    numpy_cache = None
    if args.numpy_cache is not None:
        numpy_cache: str = os.path.abspath(args.numpy_cache)
        if not numpy_cache.endswith('.npy'):
            raise argparse.ArgumentError(args.numpy_cache, message='Provided cache file does not end in .npy')
    
    #TUNING VARIABLES
    if not (0 <= args.a1min <= 1): 
        raise Exception("a1 is not in range 0 to 1")
    
    if args.min > args.max: 
        raise Exception("min nr. of partitions is higher than max")
    
    if (0 > args.chunksize):
        raise Exception('Chunksize cannot be negative')
    
    if 0 >= args.stepsize:
        raise Exception('stepsize must be a positive number larger than 0')
    
    if (args.rand == 0) is False:
        print('Setting seed to: ' + str(args.rand))
        random.seed(args.rand)
        np.random.seed(args.rand)
    
    if len(args.LList) == 0:
        raise Exception("Common contig length list has length 0! Change it ;)")

    print('Partition generation algorithm: ' + args.method)
    
    try:
        logfile = None
        if args.logfile is not None:
            logfile_path = os.path.abspath(args.logfile)
            os.remove(logfile_path)
            logfile = open(logfile_path)
            
        run(args.a1min, args.min, args.max, args.minSize, args.stepsize, args.method,\
            fasta_path, abundance_path, SCG_path, MS_path, outfile, args.numpy_cache, args.chunksize, args.LList, logfile, partition_outdir)
        
    finally:
        if logfile is not None:
            logfile.close()
    

if __name__ == "__main__":
    main()