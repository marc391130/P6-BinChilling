from io import TextIOWrapper
from EnsemblerTools import BinLogger, MergeRegulator, print_result, target_bin_3_4th_count_estimator
from ClusterDomain import Partition, Cluster, PartitionSet
from BinReaders import ContigReader, PartitionSetReader, SCGReader, ContigFilter
import numpy as np
import os
import shutil
import argparse
import random
import Constants as const
from scipy.sparse import lil_matrix, csr_matrix
import scipy.sparse as sp
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster._kmeans import euclidean_distances, KMeans, check_random_state, row_norms, MiniBatchKMeans, stable_cumsum
from sklearn.metrics import silhouette_score
from Domain import ContigData, compute_3_4_scg_count, compute_scgs_count, compute_max_scg_count
from typing import List, Dict, TypeVar, Tuple
from tqdm import tqdm
from BinChillingEnsembler import BinChillingEnsembler, Chiller, Binner
from BinEvaluator import BinEvaluator, BinRefiner
from time import time
import sys
import Assertions as Assert
import functools

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


#This has been taken from the metabinner implementation
def silhouette(X, W, label, len_weight):
    X_colsum = np.sum(X**2, axis=1)
    X_colsum = X_colsum.reshape(len(X_colsum), 1)
    W_colsum = np.sum(W**2, axis=1)
    W_colsum = W_colsum.reshape(len(W_colsum), 1)

    Dsquare = np.tile(X_colsum, (1, W.shape[0])) + np.tile(W_colsum.T, (X.shape[0], 1)) - 2 * X.dot(W.T)
    # avoid error caused by accuracy
    Dsquare[Dsquare < 0] = 0
    D = np.sqrt(Dsquare)
    aArr = D[np.arange(D.shape[0]), label]
    D[np.arange(D.shape[0]), label] = np.inf
    bArr = np.min(D, axis=1)
    tmp = (bArr - aArr) / np.maximum(aArr, bArr)
    return np.average(tmp, weights=len_weight)

def compute_partition_range(features: List[List[float]], weigths: List[float], contigs: List[ContigData], stepsize: int, min_partitions: int or None = None, max_partitions: int or None = None) -> Tuple[int, int]:
    k0 = max(compute_3_4_scg_count(contigs), 2) #Must have at least 2 clusters
    end_k = min(len(features), (k0 + max_partitions if max_partitions is not None else k0*5 + 2*stepsize+1))
    
    data = np.array(features)
    
    best_K1, best_sil_val_1 = k0, np.NINF
    best_K2, best_sil_val_2 = k0, np.NINF
    
    min_partitions = min_partitions+1 if min_partitions is not None else 1
    print(f'Searching for best K value. (Will stop when reaching {end_k} or silhouette_score decreases)')
    for k in range(k0+1, end_k+1, stepsize):
        kmeans = KMeans(n_clusters=k, init="k-means++", n_init=10, random_state=7, max_iter=300)
        labels = kmeans.fit_predict(features, sample_weight=weigths)
        # score = silhouette_score(features, labels)
        score = silhouette(data, kmeans.cluster_centers_, labels, weigths)
        print(f"k1 value of '{k}'/{end_k} got score: {score} / {best_sil_val_1}")
        if score > best_sil_val_1:
            best_sil_val_1 = score
            best_K1 = k
        else:
            break
    
    for k in range(best_K1+2*stepsize, end_k+1, stepsize):
        kmeans = KMeans(n_clusters=k, init="k-means++", n_init=10, random_state=7, max_iter=300)
        labels = kmeans.fit_predict(features, sample_weight=weigths)
        # score = silhouette_score(features, labels)
        score = silhouette(data, kmeans.cluster_centers_, labels, weigths)
        print(f"k2 value of '{k}'/{end_k} got score: {score} / {best_sil_val_2}")
        if score > best_sil_val_2:
            best_sil_val_2 = score
            best_K2 = k
        else:
            break
    
    best_K = best_K1 if best_sil_val_1 > best_sil_val_2 else best_K2
    return (k0, best_K)
                                
                                
class FeaturesDto:
    def __init__(self, index_map: Dict[int, ContigData],\
        len_weights: np.ndarray, constraint_matrix: csr_matrix or None) -> None:
        self.index_map = index_map
        self.len_weights = len_weights
        self.constraint_matrix = constraint_matrix
    def asTuple(self) -> Tuple[Dict[int, ContigData], np.ndarray, csr_matrix or None]:
        return (self.index_map, self.len_weights, self.constraint_matrix)
                                                            
def transform_contigs_to_features(items: List[ContigData], include_constraint_matrix: bool = False) \
    -> Tuple[np.ndarray, np.ndarray, np.ndarray, FeaturesDto]:
        #returns: contig index map, features, weights, constraint matrix
    index_map, comp_matrix, cov_matrix, len_weights = {}, [], [], []
    for i in range(len(items)):
        item = items[i]
        index_map[i] = item
        comp_matrix.append( item.AsNormalizedCompositionVector() )
        cov_matrix.append( item.AsNormalizedAbundanceVector() )
        # weights.append( item.avg_abundance + 0.001 )
        len_weights.append( item.contig_length )
        
    comp_matrix = np.array(comp_matrix)
    cov_matrix = np.array(cov_matrix)
    combo_matrix = np.hstack((comp_matrix, cov_matrix))
    constraint_matrix = compute_constraints(index_map) if include_constraint_matrix else None
    return (combo_matrix, comp_matrix, cov_matrix, FeaturesDto(index_map, len_weights, constraint_matrix))


def partial_seed_init(features: np.ndarray, n_clusters: int, random_state: int, scg: str, index_map: Dict[int, ContigData]):
    sorted_contigs: List[Tuple[int, ContigData]] = sorted(index_map.items(), key=lambda x: x[1].contig_length, reverse=True)
    indexes: set[int] = set()
    for index, contig in sorted_contigs:
        if scg in contig.SCG_genes:
            indexes.add(index)
    
    while len(indexes) < n_clusters:
        j = np.random.randint(0, len(features))
        indexes.add(j)
            
    matrix = np.array([features[x] for x in indexes][0:n_clusters], dtype=np.float16)
    # np.savetxt('test.txt', list_result, fmt='%.2e', newline='\n\n')
    return matrix
    
def partial_seed_init3(features: np.ndarray, n_clusters: int, random_state, seed_idx: List[int], n_local_trials=None):
    random_state = check_random_state(random_state)
    x_squared_norms = row_norms(features, squared=True)

    n_samples, n_features = features.shape

    centers = np.empty((n_clusters, n_features), dtype=features.dtype)

    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(n_clusters))

    # Pick first center randomly

    center_id = seed_idx[0]

    if sp.issparse(features):
        centers[0] = features[center_id].toarray()
    else:
        centers[0] = features[center_id]

    # Initialize list of closest distances and calculate current potential
    closest_dist_sq = euclidean_distances(
        centers[0, np.newaxis], features, Y_norm_squared=x_squared_norms,
        squared=True)

    for c, center_id in enumerate(seed_idx[1:], 1):
        if sp.issparse(features):
            centers[c] = features[center_id].toarray()
        else:
            centers[c] = features[center_id]
            # print(c, center_id)
        closest_dist_sq = np.minimum(closest_dist_sq,
                                     euclidean_distances(
                                         centers[c, np.newaxis], features, Y_norm_squared=x_squared_norms,
                                         squared=True))
    current_pot = closest_dist_sq.sum()

    # Pick the remaining n_clusters-1 points
    for c in range(len(seed_idx), n_clusters):
        # Choose center candidates by sampling with probability proportional
        # to the squared distance to the closest existing center
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(stable_cumsum(closest_dist_sq),
                                        rand_vals)
        # XXX: numerical imprecision can result in a candidate_id out of range
        np.clip(candidate_ids, None, closest_dist_sq.size - 1,
                out=candidate_ids)

        # Compute distances to center candidates
        distance_to_candidates = euclidean_distances(
            features[candidate_ids], features, Y_norm_squared=x_squared_norms, squared=True)

        # Decide which candidate is the best
        best_candidate = None
        best_pot = None
        best_dist_sq = None
        for trial in range(n_local_trials):
            # Compute potential when including center candidate
            new_dist_sq = np.minimum(closest_dist_sq,
                                     distance_to_candidates[trial])
            new_pot = new_dist_sq.sum()

            # Store result if it is the best local trial so far
            if (best_candidate is None) or (new_pot < best_pot):
                best_candidate = candidate_ids[trial]
                best_pot = new_pot
                best_dist_sq = new_dist_sq

        # Permanently add best center candidate found in local tries
        if sp.issparse(features):
            centers[c] = features[best_candidate].toarray()
        else:
            centers[c] = features[best_candidate]
        current_pot = best_pot
        closest_dist_sq = best_dist_sq

    return centers
                



#returns labels
def run_clustering_method(method: str, n_clusters: int, scg_count: Dict[str, int], matrix: np.ndarray, data:FeaturesDto, max_iter: int):
    index_map, weights, constraints = data.asTuple()
    if method == 'Kmeans':
        return KMeans(n_clusters=n_clusters, max_iter=max_iter, n_init=10, random_state=n_clusters).fit_predict(matrix, sample_weight=weights)
    if method == 'PartialSeed':
        # scg_seed = list(scg_count.keys())[ np.random.randint(0, len(scg_count)) ]
        # return KMeans(n_clusters=n_clusters, random_state=7, n_init=1,
        #             init=functools.partial(partial_seed_init, scg=scg_seed, index_map=index_map))\
        #             .fit_predict(matrix, sample_weight=weights)
        
        scg_seed = list(scg_count.keys())[ np.random.randint(0, len(scg_count)) ]
        seed_idx = [i for i, contig in index_map.items() if scg_seed in contig.SCG_genes]
        return KMeans(n_clusters=n_clusters, random_state=7, n_init=10,
                    init=functools.partial(partial_seed_init3, seed_idx=seed_idx))\
                    .fit_predict(matrix, sample_weight=weights)
    if method == 'Hierarchical':
        return AgglomerativeClustering(n_clusters=n_clusters, connectivity=constraints, memory=CAHCE_DIR).fit_predict(matrix)
    if method == 'Random':
        return KMeans(n_clusters=n_clusters, init='random', n_init=1, max_iter=max_iter).fit_predict(np.random.rand( len(matrix), 32 ) )
    raise Exception(f'method name {method} not recognized')
    

def generate_partition(partition: Partition, matrix: np.ndarray, data: FeaturesDto,\
    method: str, n_clusters: int, scg_count: Dict[str, int] , max_iter: int = 300) -> Partition:
    
    labels = run_clustering_method(method, n_clusters, scg_count, matrix, data, max_iter)
    for i in range(len(labels)):
        partition.add(str(labels[i]), data.index_map[i])
    return partition

def generate_partition_with_matrix(gamma: PartitionSet, matrix: np.ndarray, k_range: Tuple[int, int], scg_count: Dict[str, int],\
    data: FeaturesDto, method: str, max_iter: int, logger: BinLogger, partition_outdir: str or None = None):
    
    k_min, k_max = k_range
    for k in tqdm(range(k_min, k_max)):
        partition = generate_partition(gamma.create_partition(), matrix, data, method, k, scg_count, max_iter)
        if partition_outdir is not None:
            print_result(os.path.join(partition_outdir, 'Partition_'+str(k) +'.tsv'), partition)
    pass

def run_binner(a1min: float, min_partitions_gamma: int, max_partitions_gamma: int, min_contig_len:int, stepsize:int, method: str,\
    fasta_file: str, abundance_file: str, scg_file: List[str], ms_file: List[str], output_file: str, numpy_cache: str, chunksize: int,\
        logger: BinLogger, max_iter: int = 300, partition_outdir: str or None = None):
    
    start_time = time()
    scg_reader = SCGReader(scg_file, ms_file, logger=logger)
    filter = ContigFilter(min_contig_len)
    contig_reader = ContigReader(fasta_file, scg_reader, depth_file=abundance_file,\
        numpy_file=numpy_cache, enable_analyse_contig_comp=True, logger=logger)
    
    contigs = [x for x in list(contig_reader.read_file_fast(load_SCGs=True).values()) if filter.predicate(x)] 
    if any( (len(x.SCG_genes) > 0 for x in contigs) ) is False:
        raise Exception('The set of contigs contains no SCGs')
    if any( (x.__has_analysis__() for x in contigs) ) is False:
        raise Exception('The contigs did not contain an analysis of composition. '+\
            'This might be because the cache file is generated using fast analysis. '+\
            'Try deleting cache file, or provide a different cache file path')
    
    gamma = PartitionSet(contigs)
    feature_data = transform_contigs_to_features(contigs, include_constraint_matrix= (method == 'Hierarchical'))
    combo_matrix, comp_matrix, cov_matrix, data_dto = feature_data
    
    scg_count = compute_scgs_count(contigs)
    try:
        k_min, k_max = compute_partition_range(comp_matrix, data_dto.len_weights, contigs, stepsize, min_partitions_gamma, max_partitions_gamma)
        k_max = max(k_min + min_partitions_gamma, k_max)
        k_range = (k_min, k_max)
        
        logger.log(f'Generating partitions with {k_min} to {k_max} bins using features matricies...')
        # logger.log(f'Generating partitions using combo matrix...')
        # generate_partition_with_matrix(gamma, combo_matrix, k_range, scg_count,\
        #     data_dto, method, max_iter, logger, partition_outdir)
        
        logger.log(f'Generating partitions using composition_matrix...')
        generate_partition_with_matrix(gamma, comp_matrix, k_range, scg_count,\
            data_dto, method, max_iter, logger, partition_outdir)
        
        logger.log(f'Generating partitions using composition_matrix...')
        generate_partition_with_matrix(gamma, comp_matrix, k_range, scg_count,\
            data_dto, method, max_iter, logger, partition_outdir)
        
        logger.log(f'Generating partitions using composition_matrix...')
        generate_partition_with_matrix(gamma, comp_matrix, k_range, scg_count,\
            data_dto, method, max_iter, logger, partition_outdir)
        
        # logger.log(f'Generating partitions using abundance matrix...')
        # generate_partition_with_matrix(gamma, cov_matrix, k_range, scg_count,\
        #     data_dto, method, max_iter, logger, partition_outdir)
        
        # generate_partition_with_matrix(gamma, combo_matrix, k_range, scg_count,\
        #     data_dto, 'Kmeans', max_iter, logger, partition_outdir)
        
        # generate_partition_with_matrix(gamma, comp_matrix, k_range, scg_count,\
        #     data_dto, 'Kmeans', max_iter, logger, partition_outdir)
        
        # generate_partition_with_matrix(gamma, cov_matrix, k_range, scg_count,\
        #     data_dto, 'Kmeans', max_iter, logger, partition_outdir)
        #end_loop
    finally:
        shutil.rmtree(os.path.abspath(CAHCE_DIR), ignore_errors=True)
        
    
    bin_evaluator = BinEvaluator(scg_reader.read_MS_scgs())
    bin_refiner = BinRefiner(bin_evaluator, (1.0 / len(gamma)), logger)
    chiller = Chiller(a1min, 1.0, MergeRegulator(a1min), 0.02, logger)
    binner = Binner(bin_refiner, bin_evaluator, logger=logger)
    ensembler = BinChillingEnsembler(chiller, binner, bin_evaluator, chunksize=chunksize, target_clusters_est=target_bin_3_4th_count_estimator, logger=logger)

    final_partition = ensembler.ensemble(gamma)
    print_result(output_file, final_partition)
    
    logger.log(f"Completed binning in time: {(time() - start_time):0.02f}s")
