from BinEvaluator import BinEvaluator
from Cluster import Cluster, Partition, PartitionSet
import numpy as np
from typing import List, Dict, Tuple, Callable
from tqdm import tqdm
from ClusterSimilarityMatrix import SparseClustserSimularity
from ClusterSimilarityMatrix import cluster_simularity
from multiprocessing import cpu_count, Pool
from Domain import ContigData
from MemberSimularityMatrix import MemberSimularityMatrix

def print_result(file_path: str, parititon: Partition[ContigData]):
    with open(file_path, 'w') as file:
        cluster_lst = list(parititon.values())
        for cluster_idx in range(len(cluster_lst)):
            for item in cluster_lst[cluster_idx]:
                file.write(f"{cluster_idx+1}\t{item.name}\n")


class MergeRegulator:
    def __init__(self, a1_min: float) -> None:
        self.a1_min = a1_min
        self.target_clusters = 0
        self.result = []
    
    def set_context(self, gamma: PartitionSet, target_clusters: int):
        self.target_clusters = target_clusters
        self.__context__ = gamma
    
    def evaluate(self, alpha1: float, cluster_matrix: SparseClustserSimularity, merged_clusters: List[Cluster]) \
        -> Tuple[bool, List[Cluster]]:
        if alpha1 < self.a1_min: return True, self.result
        self.result = cluster_matrix.get_available_clusters() + merged_clusters
        
        return len(self.result) < self.target_clusters, self.result
    
    def get_merge_result(self) -> List[Cluster]:
        self.target_clusters = 0
        self.__context__ = None
        result = self.result
        self.result = []
        return result


        
def handle_estimate_target_clusters(gamma: PartitionSet,  taget_clusters_est: int or Callable[[PartitionSet], int]) -> int:
        if isinstance(taget_clusters_est, int):
            return taget_clusters_est
        elif callable(taget_clusters_est):
            return int(taget_clusters_est(gamma))
        else:
            return int(max([len(partition) for partition in gamma]))

def target_bin_3_4th_count_estimator(gamma: PartitionSet) -> int:
    partition_ln = [len(partition) for partition in gamma]
    average = sum(partition_ln) / len(partition_ln)
    third = (max(partition_ln) - average ) / 2
    return int(average + third)

def sort_merged_cluster_singlethread(cluster_matrix: SparseClustserSimularity, merged_lst: List[Cluster]) -> float:
    max_simularity = -1
    for merged_cluster in merged_lst:
        cluster_matrix.add_cluster(merged_cluster)
    
    for merged_cluster in tqdm(merged_lst):
        for cluster in cluster_matrix.__all_clusters__:
            if merged_cluster is cluster: continue
            similarity = cluster_simularity(merged_cluster, cluster, cluster_matrix.total_item_count)
            if similarity > cluster_matrix.min_value:
                cluster_matrix.set_value(merged_cluster, cluster, similarity)
            #elif similarity < alpha1: 
            max_simularity = max(max_simularity, similarity)
    return max_simularity

def sort_merged_cluster_multithread(cluster_matrix: SparseClustserSimularity, merged_lst: List[Cluster],\
    threads: int = cpu_count(), chunksize:int = 1) -> float:
    max_simularity = -1
    for c in merged_lst: cluster_matrix.add_cluster(c)
    index_lst = list(cluster_matrix.__all_clusters__)
    parameters = [(cluster, cluster_matrix.total_item_count, cluster_matrix.min_value) for cluster in merged_lst]
    global merged_cluster_presort_cache
    
    try:
        merged_cluster_presort_cache = index_lst
        
        with Pool(threads) as p:
            for index, res, max_sim in tqdm(p.imap_unordered(partial_sort_merge_cache, parameters, chunksize=chunksize), total=len(parameters)):
                max_simularity = max(max_simularity, max_sim)
                for index2, sim in res:
                    cluster_matrix.set_value(index_lst[index], index_lst[index2], sim)
    finally:
        merged_cluster_presort_cache = None
        
    
    return max_simularity

def partial_sort_merge(tup: Tuple[Cluster, List[Cluster], int, float]) -> Tuple[int, List[Tuple[int, float]], float]:
    merged_cluster, all_clusters, total_item_count, min_value = tup
    own_index, max_similarity, result = -1, -1, []
    for cluster_idx in range(len(all_clusters)):
        cluster = all_clusters[cluster_idx]
        if merged_cluster is cluster:
            own_index = cluster_idx           
            continue
        similarity = cluster_simularity(merged_cluster, cluster, total_item_count)
        max_similarity = max(similarity, max_similarity)
        if similarity > min_value:
            result.append( (cluster_idx, similarity) )
    return own_index, result, max_similarity


merged_cluster_presort_cache: List[Cluster] = None
def partial_sort_merge_cache(tup: Tuple[Cluster, int, float]) -> Tuple[int, List[Tuple[int, float]], float]:
    merged_cluster, total_item_count, min_value = tup
    global merged_cluster_presort_cache
    all_clusters = merged_cluster_presort_cache
    own_index, max_similarity, result = -1, -1, []
    for cluster_idx in range(len(all_clusters)):
        cluster = all_clusters[cluster_idx]
        if merged_cluster is cluster:
            own_index = cluster_idx           
            continue
        similarity = cluster_simularity(merged_cluster, cluster, total_item_count)
        max_similarity = max(similarity, max_similarity)
        if similarity > min_value:
            result.append( (cluster_idx, similarity) )
    return own_index, result, max_similarity


def MergeClusters(alpha1: float, clusters: List[Cluster], total_elements)\
    -> Tuple[List[Cluster], float]: #tuple of (all available cluster, next_higext_similarity) 
    result_clusters = []
    child_merged_set = set() #has been merged, aka skip if in this
    max_similarity = -1

    for i in tqdm(range(len(clusters))):
        cluster1 = clusters[i]
        is_merged = False
        if cluster1 in child_merged_set:
            continue

        for j in range(i+1, len(clusters)):
            cluster2 = clusters[j]

            if cluster2 in child_merged_set or cluster1.SamePartitionAs(cluster2):
                continue
                
            # if cluster_sim_matrix[cluster1, cluster2] >= alpha1:
            similarity = cluster_simularity(cluster1, cluster2, total_elements)
            if similarity >= alpha1:
                #merge the two clusters
                new_cluster = Cluster.merge(cluster1, cluster2)
                
                #add them to the skip set
                child_merged_set.add(cluster1)
                child_merged_set.add(cluster2)
                
                #append to results
                result_clusters.append(new_cluster)
                
                is_merged = True
                break
            else:
                max_similarity = max(max_similarity, similarity)
            
        if is_merged is False:
            result_clusters.append(cluster1)
            
        
    return result_clusters, max_similarity



#returns (item_id, s_x, cluster_index)
def __partial_cluster_certainty_degree__(\
    tuple : Tuple[int, np.matrix] ) -> Tuple[int, float, int]:
    # item, id, cluster_dct, ensembler = tuple
    item_index, matrix = tuple
    
    cluster_index = np.argmax(matrix[item_index])
    similarity = matrix[item_index, cluster_index]
    return (item_index, similarity, cluster_index)