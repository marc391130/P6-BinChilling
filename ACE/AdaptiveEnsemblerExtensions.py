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

class QualityMeasuerer:
    def calculate_quality(self, cluster: Cluster, partition_count: int, simularity_matrix: MemberSimularityMatrix) -> float:
        if len(cluster) == 0:
            return 0.0
        
        cluster_sim_dct = simularity_matrix.get_column(cluster)
        mean = sum(cluster_sim_dct.values()) / len(cluster_sim_dct)
        
        sum_value = sum([ (sim - mean)**2 for sim in cluster_sim_dct.values() ])
        
        return sum_value / len(cluster)
    
    def calculate_speculative_quality(self, initial_quality: float, include_item: object, \
        cluster: Cluster, gamma: PartitionSet, simularity_matrix: MemberSimularityMatrix) -> float:
        if include_item in cluster:
            return initial_quality

        mean = simularity_matrix.cluster_mean(cluster)
        total_var = initial_quality * len(cluster)
        simularity = (simularity_matrix.getEntry(include_item, cluster) - mean)**2
        

        return (total_var + simularity) / (len(cluster) +1)
    
    def calculate_excluding_quality(self, cluster: Cluster, exclude_item: object, partition_count: int, similarity_maxtrix: MemberSimularityMatrix) -> float:

        if len(cluster) == 0:
            return 0

        if len(cluster) == 1 and exclude_item in cluster:
            return 0

        if exclude_item not in cluster:
            return self.calculate_quality(cluster, partition_count, similarity_maxtrix)

        all_member_similarity = similarity_maxtrix.get_column(cluster)

        fuckery = all_member_similarity.get(exclude_item)

        all_member_similarity.pop(exclude_item, None)

        mean = sum(all_member_similarity.values()) / (len(cluster) - 1)

        fuckery = (fuckery - mean)**2

        quality = sum([(value - mean)**2 for value in all_member_similarity.values()])
        return (quality + fuckery / len(cluster)) - (quality / (len(cluster) - 1))

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

class AssignRegulator:
    def __init__(self, quality_measure: QualityMeasuerer, logger: Callable[[str], None] = None) -> None:
        self.logger = logger
        self.quality_measure = quality_measure

    def log(self, string: str) -> None:
        if self.logger is not None:
            self.logger(string)

    def assign_items(self, candidate_clusters: List[Cluster], totally_certain_lst: List[object], certain_lst: List[object], \
        uncertain_lst: List[object], totally_uncertain_map: Dict[object, object], gamma: PartitionSet,\
            similarity_matrix:MemberSimularityMatrix, lost_items: List[object]) -> List[Cluster]:

        self.log("Assigning totally certain objects...")
        candidate_clusters = self.__assign_certains__(totally_certain_lst, candidate_clusters)
        
        self.log("Assign certain objects...")
        candidate_clusters = self.__assign_certains__(certain_lst, candidate_clusters)
        
        self.log("Assign uncertain objects")
        candidate_clusters = self.__assign_uncertains__(uncertain_lst, candidate_clusters, gamma, similarity_matrix)
    
        candidate_clusters = self.__assign_lost_objects__(candidate_clusters, totally_uncertain_map, lost_items)
        
        return candidate_clusters


    def __assign_certains__(self, certain_lst: List[Tuple[object, Cluster]],\
        candidate_clusters: List[Cluster]) -> List[Cluster]:
        
        for item_cluster in tqdm(certain_lst):
            #done so type hinting can actually be done. Damn you tqdm
            item: object = item_cluster[0]
            cluster: Cluster = item_cluster[1]
            
            
            #remove item from all other clusters in candidate clusters
            for can_cluster in candidate_clusters:
                can_cluster.remove(item)
            #add it back into best cluster
            cluster.add(item)
        return candidate_clusters

    def __assign_uncertains__(self, uncertain_item_lst: List, candidate_clusters: List[Cluster],\
        gamma: PartitionSet, simularity_matrix: MemberSimularityMatrix) -> List[Cluster]:
        
        self.log("Calculate initial quality...")
        initial_quality = {}
        for cluster in tqdm(candidate_clusters):
            quality = self.quality_measure.calculate_quality(cluster, len(gamma), simularity_matrix)
            initial_quality[cluster] = quality
        
        self.log("Assigning uncertain clusters...")
        for item in tqdm(uncertain_item_lst):
            #this is cursed
            #takes all candidate clusters, makes a tuple of (cluster, speculative_quality)
            #then finds the cluster and new quality with the smallest delta quality,
            # by taking the abselute value of initial quality - speculative quality
            min_changed_cluster, new_quality = min([(cluster, self.quality_measure.calculate_speculative_quality( \
                initial_quality[cluster], item, cluster, gamma, simularity_matrix)) \
                for cluster in candidate_clusters], key= lambda x: abs(initial_quality[x[0]] - x[1]) )
            initial_quality[min_changed_cluster] = new_quality
            
            for cluster in candidate_clusters:
                cluster.remove(item)
                
            min_changed_cluster.add(item)
        return candidate_clusters

    def __assign_lost_objects__(self, candidate_clusters: List[Cluster], assosiation_map: Dict[object, object],\
        lost_items: List[object]) -> List[Cluster]:
        if len(assosiation_map) > 0:
            self.log('Trying to assigning remaining totally uncertain objects using co-assosiation..')
            for item, assosiate_item in tqdm(assosiation_map.items()):
                found = False
                for cluster in candidate_clusters:
                    if assosiate_item in cluster:
                        cluster.append(item)
                        found = True
                        break
                if not found: #This eliminates the circular reference problem 
                    isolated_cluster = Cluster()
                    isolated_cluster.add(item)
                    candidate_clusters.append(isolated_cluster)
                    
        if len(lost_items) > 0:
            self.log(f"Adding remaining '{len(lost_items)}' items to isolated clusters...")
            for item in lost_items:
                isolated_cluster = Cluster()
                isolated_cluster.add(item)
                candidate_clusters.append(isolated_cluster)
        return candidate_clusters

        
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
    parameters = [(cluster, index_lst, cluster_matrix.total_item_count, cluster_matrix.min_value) for cluster in merged_lst]
    
    with Pool(threads) as p:
        for index, res, max_sim in tqdm(p.imap_unordered(partial_sort_merge, parameters, chunksize=chunksize), total=len(parameters)):
            max_simularity = max(max_simularity, max_sim)
            for index2, sim in res:
                cluster_matrix.set_value(index_lst[index], index_lst[index2], sim)
    
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