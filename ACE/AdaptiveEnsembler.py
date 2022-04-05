from multiprocessing import cpu_count, Pool
from typing import Callable, Dict, List, Tuple

from more_itertools import partition
from Cluster import Cluster, Contig, Partition, PartitionSet
import numpy as np
from tqdm import tqdm
from itertools import islice
import Assertions as Assert
from sys import maxsize as MAXSIZE
from time import time
from MemberSimularityMatrix import MemberMatrix, MemberSimularityMatrix
from Logger import Logger
from ClusterSimilarityMatrix import ClusterSimilarityMatrix
from io import TextIOWrapper

THREAD_COUNT = min(cpu_count(), 8) 

MERGED_CLUSTER = -1
TOTALLY_CERTAIN, CERTAIN, UNCERTAIN, TOTALLY_UNCERTAIN = 3, 2, 1, 0

class Ensembler:
    def ensemble(self, data: PartitionSet) -> Partition:
        pass
    
class QualityMeasuerer:
    def calculate_quality(self, cluster: Cluster, member_matrix: MemberMatrix,\
        simularity_matrix: MemberSimularityMatrix) -> float:
        if len(cluster) == 0:
            return 0.0
        
        mean = simularity_matrix.Cluster_Mean(cluster)
        sum_value = sum([ pow(simularity_matrix.GetEntry(cluster, item) - mean, 2) for item in cluster])
        
        return sum_value / len(cluster)
    
    def calculate_speculative_quality(self, initial_quality: float, item: object, cluster: Cluster,\
        member_matrix: MemberMatrix, simularity_matrix: MemberSimularityMatrix) -> float:
        if item in cluster:
            return initial_quality
        
        mean = simularity_matrix.Cluster_Mean(cluster)
        member_val = member_matrix.getEntry(cluster, item)
        div = simularity_matrix.Item_sum(item)
        new_simularity = (member_val +1) / (div+1)
        value = pow(new_simularity - mean, 2)
        
        total_val = initial_quality * len(cluster)
        
        return (total_val + value) / (len(cluster) +1)

def target_bin_3_4th_count_estimator(gamma: PartitionSet) -> int:
    cluster_ln = [len(partition) for partition in gamma]
    average = sum(cluster_ln) / len(cluster_ln)
    third = (max(cluster_ln) - average ) / 2
    return int(average + third)
    

class AdaptiveClusterEnsembler(Ensembler):
    def __init__(self, initial_alpha1_thredshold: float = 0.8, initial_delta_aplha: float = 0.1,\
        alpha1_min: float = 0.6, alpha2: float = 0.2, taget_clusters_est: int or Callable[[PartitionSet], int] = None, quality_measurer: QualityMeasuerer = None, logfile: TextIOWrapper = None):
        self.alpha1_thredshold = initial_alpha1_thredshold
        self.delta_alpha = initial_delta_aplha
        self.aplha1_min = alpha1_min
        self.alpha2 = alpha2
        self.logfile = logfile
        self.quality_measure = quality_measurer if quality_measurer is not None else QualityMeasuerer()
        self.start_time: float = None
        self.taget_clusters_est = taget_clusters_est
    
    
    def __target_clusters__(self, gamma: PartitionSet) -> int:
        if isinstance(self.taget_clusters_est, int):
            return self.taget_clusters_est
        elif callable(self.taget_clusters_est):
            return int(self.taget_clusters_est(gamma))
        else:
            return int(max([len(partition) for partition in gamma]))
    
    def log(self, string: str) -> None:
        print(string, file=self.logfile)

            
    def ensemble(self, gamma: PartitionSet) -> Partition:
        
        self.start_time = time()
        
        all_items = list(gamma.get_all_elements().keys())
        all_clusters = gamma.get_all_clusters()
        
        self.log("Building Cluster simularity matrix...")
        clusterMatrix = ClusterSimilarityMatrix.Build(all_clusters, len(all_items))
        target_clusters = self.__target_clusters__(gamma)
        
        memberMatrix: MemberMatrix = None
        
        self.log("started 2.1")
        lambda_len, available_clusters = 0, all_clusters
        while True:
            available_clusters = MergeClusters(self.alpha1_thredshold, clusterMatrix, all_clusters)
            lambda_len = len(available_clusters)
            self.log(f"Testing alhpa1 value of {self.alpha1_thredshold}, found {lambda_len} / { target_clusters} ...")
            if lambda_len >= target_clusters:
                # memberMatrix = MemberMatrix(available_clusters, all_items)
                break
            else:
                self.alpha1_thredshold += self.delta_alpha
        
        
        self.log("started 2.2")
        i, start_22 = 0, time()
        def log_22() -> None:
            self.log(f"iteration {i:05d}, alpha1: {self.alpha1_thredshold:.4f} / { self.aplha1_min }, clusters: {lambda_len} / { target_clusters }, Time: {(time() - start_22):0.02f}s")
        log_22()
        while lambda_len >= target_clusters:
            i += 1
            clusterMatrix = ClusterSimilarityMatrix.BuildFrom(available_clusters, clusterMatrix)
            self.alpha1_thredshold = clusterMatrix.max_simularity_np()
            log_22()
            if self.alpha1_thredshold < self.aplha1_min:
                break
            else:
                new_available_clusters = MergeClusters(self.alpha1_thredshold, clusterMatrix, available_clusters)
        
            if len(new_available_clusters) < target_clusters:
                break
            else:
                available_clusters = new_available_clusters
                lambda_len = len(available_clusters)
        
        self.log("started step 2.3, building simularity matrix...")
        memberMatrix = MemberMatrix.build(available_clusters, all_items)
        similarity_matrix = memberMatrix.BuildSimularityMatrix(available_clusters)
        certain_clusters = [cluster for cluster in tqdm(available_clusters) if similarity_matrix.Cluster_Max(cluster) > self.alpha2]
        
        self.log(f"Found {len(certain_clusters)} clusters with certain objects")
        candidate_clusters, non_candidate_clusters = None, None
        if len(certain_clusters) == target_clusters:
            candidate_clusters = certain_clusters
            non_candidate_clusters = [cluster for cluster in available_clusters if cluster not in candidate_clusters]
        else: 
            cluster_certainty_lst = [ (cluster, similarity_matrix.Cluster_Mean(cluster)) for cluster in available_clusters ]
            cluster_certainty_lst = list(sorted(cluster_certainty_lst, key = lambda item: item[1], reverse=True))
            candidate_clusters: List[Cluster] = [x[0] for x in islice(cluster_certainty_lst, target_clusters)]
            non_candidate_clusters = [x[0] for x in islice(cluster_certainty_lst, target_clusters, None)]
            self.alpha2 = similarity_matrix.Cluster_Max(candidate_clusters[len(candidate_clusters)-1])
        
        #calculate new membership and simularity based on candidate and non-candidate
        candidate_memberMatrix = MemberMatrix.build(candidate_clusters, all_items)
        non_candidate_memberMatrix = MemberMatrix.build(non_candidate_clusters, all_items)
        similarity_matrix = MemberSimularityMatrix.RefineTo(similarity_matrix, candidate_clusters)
        
        partition = self.assign_item_to_one_cluster(gamma,\
            candidate_clusters, non_candidate_clusters, similarity_matrix, candidate_memberMatrix, non_candidate_memberMatrix)
        
        self.log(f"Finished in time {(time() - self.start_time):0.02f}s")
        self.start_time = None
        
        return partition
    
    
    def assign_item_to_one_cluster(self, gamma: PartitionSet, candidate_clusters: List[Cluster], non_candidate_clusters: List[Cluster],\
        similarity_matrix: MemberSimularityMatrix, membership_matrix: MemberMatrix, non_candidate_membership_matrix) -> Partition:
        
        all_items = gamma.get_all_elements()
        
        self.log("Determining object certainty...")
        totally_certain_lst, certain_lst, uncertain_lst, totally_uncertain_lst \
            = self.identify_object_certainy(all_items, similarity_matrix)
        
        self.log("\nFound: ")
        self.log(f"Totally certain objects: {len(totally_certain_lst)}")
        self.log(f"Certain objects: {len(certain_lst)}")
        self.log(f"Uncertain objects: {len(uncertain_lst)}")
        self.log(f"Totally uncertain objects: {len(totally_uncertain_lst)}")
        self.log('\n')
        
        
        
        if len(totally_uncertain_lst) > 0:
            self.log("reidentifying totally uncertain items...")
            totally_certain_lst_2, certain_lst_2, uncertain_lst_2 = self.reidentify_totally_uncertain_item(totally_uncertain_lst, candidate_clusters,\
                non_candidate_membership_matrix, similarity_matrix, gamma)
            
            
            totally_certain_lst += totally_certain_lst_2
            certain_lst += certain_lst_2
            uncertain_lst += uncertain_lst_2
            
            self.log("\nObject certainty after totally uncertain items have been reidentified: ")
            self.log(f"Totally certain objects: {len(totally_certain_lst)}")
            self.log(f"Certain objects: {len(certain_lst)}")
            self.log(f"Uncertain objects: {len(uncertain_lst)}")
            self.log('\n')
            
        else:
            self.log("no totally uncertain objects found, skipping reidentification step")
        
        
        self.log("Assigning totally certain objects...")
        candidate_clusters = self.assign_certains_objects(totally_certain_lst, candidate_clusters)
        
        self.log("Assign certain objects...")
        candidate_clusters = self.assign_certains_objects(certain_lst, candidate_clusters)
        
        self.log("Assign uncertain objects")
        candidate_clusters = self.assign_uncertain_objects(uncertain_lst, candidate_clusters, similarity_matrix, membership_matrix)
        
                        
        return self.build_final_partition(gamma, candidate_clusters)
    
    
    def reidentify_totally_uncertain_item(self, totally_uncertain_item_lst: List, candidate_clusters: List[Cluster],\
        non_can_membership: MemberMatrix, simularity_matrix: MemberSimularityMatrix, gamma: PartitionSet) -> Tuple[List, List, List]:
        
        coassosiation_matrix = non_can_membership.buildCoAssosiationMatrix(gamma)
        def recalculate_simularity(item: object) -> None:
            Assert.assert_item_in_list(totally_uncertain_item_lst, item)
            for cluster in candidate_clusters:
                v = non_can_membership.average_common_neighbors(coassosiation_matrix, item, cluster)
                simularity_matrix[item, cluster] = v
            return None
        
        for totally_uncertain_item in tqdm(totally_uncertain_item_lst):
            recalculate_simularity(totally_uncertain_item)
        
        self.log("Reclassifying totally uncertain objects...")
        tot_certain, certain, uncertain, tot_uncertain = self.identify_object_certainy(totally_uncertain_item_lst, simularity_matrix)
        Assert.assert_equal(len(tot_uncertain), 0)
        
        return tot_certain, certain, uncertain
        
        
    
    def assign_uncertain_objects(self, uncertain_item_lst: List, candidate_clusters: List[Cluster],\
        simularity_matrix: MemberSimularityMatrix, membership_matrix: MemberMatrix) -> List[Cluster]:
        
        self.log("Calculate initial quality...")
        initial_quality = {}
        for cluster in tqdm(candidate_clusters):
            quality = self.quality_measure.calculate_quality(cluster, membership_matrix, simularity_matrix)
            initial_quality[cluster] = quality
        
        self.log("Assigning uncertain clusters...")
        for item in tqdm(uncertain_item_lst):
            #this is cursed
            #takes all candidate clusters, makes a tuple of (cluster, speculative_quality)
            #then finds the cluster and new quality with the smallest delta quality,
            # by taking the abselute value of initial quality - speculative quality
            min_changed_cluster, new_quality = min([(cluster, self.quality_measure.calculate_speculative_quality( \
                initial_quality[cluster], item, cluster, membership_matrix, simularity_matrix)) \
                for cluster in candidate_clusters], key= lambda x: abs(initial_quality[x[0]] - x[1]) )
            initial_quality[min_changed_cluster] = new_quality
            
            for cluster in candidate_clusters:
                cluster.remove(item)
                
            min_changed_cluster.append(item)
        return candidate_clusters
    
    def assign_certains_objects(self, certain_lst: List[Tuple[object, Cluster]],\
        candidate_clusters: List[Cluster]) -> List[Cluster]:
        
        for item_cluster in tqdm(certain_lst):
            #done so type hinting can actually be done. Damn you tqdm
            item: object = item_cluster[0]
            cluster: Cluster = item_cluster[1]
            
            
            #remove item from all other clusters in candidate clusters
            for can_cluster in candidate_clusters:
                can_cluster.remove(item)
            #add it back into best cluster
            cluster.append(item)
        return candidate_clusters
    
        
            
    def identify_object_certainy(self, items: List, similarity_matrix: MemberSimularityMatrix) \
        -> Tuple[List[Tuple[object, Cluster]], List[Tuple[object, Cluster]], List[object], List[object]]:
        totally_certain_lst, certain_lst, uncertain_lst, totally_uncertain_lst  = [], [], [], []
    
        # np.savetxt("matrix.txt", similarity_matrix.matrix, fmt='%f', newline='\n\n')
        # np.savetxt("matrix_sum_0.txt", similarity_matrix.matrix.sum(axis=0), delimiter=',', fmt='%f', newline='\n\n')
        # np.savetxt("matrix_sum_1.txt", similarity_matrix.matrix.sum(axis=1), delimiter=',', fmt='%f', newline='\n\n')
        # np.savetxt("matrix_max_0.txt", similarity_matrix.matrix.max(axis=0), delimiter=',', fmt='%f', newline='\n\n')
        # np.savetxt("matrix_max_1.txt", similarity_matrix.matrix.max(axis=1), delimiter=',', fmt='%f', newline='\n\n')
        reverse_cluster_index_map = {value: key for key, value in similarity_matrix.cluster_index_map.items()}
        reverse_item_index_map = {value: key for key, value in similarity_matrix.item_index_map.items()}

        parameters = []
        for item in items:
            parameters.append((similarity_matrix.item_index_map[item], similarity_matrix.matrix))
        
        with Pool(THREAD_COUNT) as p:
            r = list(tqdm(p.imap(__partial_cluster_certainty_degree__, parameters), total=len(parameters)))
            for item_id, similarity, cluster_index in r:
                item_cluster_tuple = (reverse_item_index_map[item_id], reverse_cluster_index_map[cluster_index])
                
                if similarity >= 1:
                    totally_certain_lst.append(item_cluster_tuple)
                elif similarity > self.alpha2 and similarity < 1:
                    certain_lst.append(item_cluster_tuple)
                elif similarity <= self.alpha2 and similarity > 0:
                    uncertain_lst.append(item_cluster_tuple[0]) #only add item
                elif similarity == 0:
                    totally_uncertain_lst.append(item_cluster_tuple[0]) #only add item
                else:
                    raise Exception("something went wrong")
        return totally_certain_lst, certain_lst, uncertain_lst, totally_uncertain_lst
        
    def build_final_partition(self, gamma: PartitionSet, candidate_clusters: List[Cluster]):
        partition = Partition(list(gamma.get_all_elements().keys()))
        
        self.log("Building final partition from candidate clusters")
        for cluster in tqdm(candidate_clusters):
            # Assert.assert_list_nonempty(cluster)
            for item in cluster:
                partition.add(str(cluster), item)
                
        self.log(f"Found {len(partition)} total clusters")
        return partition
            
        
        
def MergeClusters(alpha1: float, cluster_sim_matrix: ClusterSimilarityMatrix, clusters: List[Cluster])\
    -> Tuple[List[Cluster], set]: #tuple of (all available cluster, set of merged clusters) 
    result_clusters = []
    child_merged_set = set() #has been merged, aka skip if in this


    for i in range(len(clusters)):
        cluster1 = clusters[i]
        is_merged = False
        if cluster1 in child_merged_set:
            continue

        for j in range(i, len(clusters)):
            cluster2 = clusters[j]

            if cluster1 is cluster2 or cluster1.SamePartitionAs(cluster2) or cluster2 in child_merged_set:
                continue
                
            if cluster_sim_matrix[cluster1, cluster2] >= alpha1:
                #merge the two clusters
                new_cluster = Cluster.merge(cluster1, cluster2)
                
                #add them to the skip set
                child_merged_set.add(cluster1)
                child_merged_set.add(cluster2)
                
                #append to results
                result_clusters.append(new_cluster)
                
                is_merged = True
                break
            
        if is_merged is False:
            result_clusters.append(cluster1)
            
        
    return result_clusters



#returns (item_id, s_x,cluster_index)
def __partial_cluster_certainty_degree__(\
    tuple : Tuple[int, np.matrix] ) -> Tuple[int, float, int]:
    # item, id, cluster_dct, ensembler = tuple
    item_index, matrix = tuple
    
    cluster_index = np.argmax(matrix[item_index])
    similarity = matrix[item_index, cluster_index]
    return (item_index, similarity, cluster_index)


def __assign_to_certainty__(item: object, simularity: float, alpha2: int,\
    totally_certain: List, certain: List, uncertain: List, totally_uncertain: List) -> None:
    
    if simularity >= 1:
        totally_certain.append(item)
    elif alpha2 < simularity < 1:
        certain.append(item)
    elif simularity <= alpha2 and simularity != 0:
        uncertain.append(item)
    elif simularity == 0:
        totally_uncertain.append(item)
    else:
        raise Exception("nay nayed")