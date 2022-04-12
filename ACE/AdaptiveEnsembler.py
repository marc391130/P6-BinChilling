from functools import partialmethod
from multiprocessing import cpu_count, Pool
from typing import Callable, Dict, List, Tuple
from Cluster import Cluster, Partition, PartitionSet
import numpy as np
from tqdm import tqdm
from itertools import islice
import Assertions as Assert
from sys import maxsize as MAXSIZE
from time import time
from MemberSimularityMatrix import CoAssosiationMatrix, MemberMatrix, MemberSimularityMatrix
from ClusterSimilarityMatrix import ClusterSimilarityMatrix
from io import TextIOWrapper

__global_disable_tqdm = False
tqdm.__init__ = partialmethod(tqdm.__init__, disable=__global_disable_tqdm)

THREAD_COUNT = min(cpu_count(), 8) 

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
    partition_ln = [len(partition) for partition in gamma]
    average = sum(partition_ln) / len(partition_ln)
    third = (max(partition_ln) - average ) / 2
    return int(average + third)
    

class AdaptiveClusterEnsembler(Ensembler):
    def __init__(self, initial_alpha1_thredshold: float = 0.8, initial_delta_aplha: float = 0.1,\
        alpha1_min: float = 0.6, alpha2: float = 0.2, \
        taget_clusters_est: int or Callable[[PartitionSet], int] = None, quality_measurer: QualityMeasuerer = None,\
        logfile: TextIOWrapper = None, should_log: bool = True, threads: int = None):
        self.alpha1_thredshold = initial_alpha1_thredshold
        self.delta_alpha = initial_delta_aplha
        self.aplha1_min = alpha1_min
        self.alpha2 = alpha2
        self.logfile = logfile
        self.quality_measure = quality_measurer if quality_measurer is not None else QualityMeasuerer()
        self.taget_clusters_est = taget_clusters_est
        self.should_log = should_log #Should log to console
        self.thread_count = min(threads, THREAD_COUNT) if threads is not None else THREAD_COUNT
    
    
    def __calc_target_clusters__(self, gamma: PartitionSet) -> int:
        if isinstance(self.taget_clusters_est, int):
            return self.taget_clusters_est
        elif callable(self.taget_clusters_est):
            return int(self.taget_clusters_est(gamma))
        else:
            return int(max([len(partition) for partition in gamma]))
    
    def log(self, string: str) -> None:
        if self.should_log:
            print(string)
        if self.logfile is not None:
            print(string, file=self.logfile)

            
    def ensemble(self, gamma: PartitionSet) -> Partition:
        
        start_time = time()
        alpha1, alpha2 = self.alpha1_thredshold, self.alpha2
        
        all_items = list(gamma.get_all_elements().keys())
        all_clusters = gamma.get_all_clusters()
        
        self.log("Building Cluster simularity matrix...")
        clusterMatrix = ClusterSimilarityMatrix.Build(all_clusters, len(all_items))
        target_clusters = self.__calc_target_clusters__(gamma)
        
        memberMatrix: MemberMatrix = None
        
        self.log("Merging initial clusters (step 2.1)")
        lambda_len, available_clusters = 0, all_clusters
        while True:
            available_clusters = MergeClusters(alpha1, clusterMatrix, all_clusters)
            lambda_len = len(available_clusters)
            self.log(f"Testing alhpa1 value of {alpha1}, found {lambda_len} / { target_clusters} ...")
            if lambda_len >= target_clusters:
                # memberMatrix = MemberMatrix(available_clusters, all_items)
                break
            else:
                alpha1 += self.delta_alpha
        
        
        self.log("Merging similar clusters (step 2.2)")
        i, start_22 = 0, time()
        def log_22() -> None:
            self.log(f"iteration {i:05d}, alpha1: {alpha1:.4f} / { self.aplha1_min }, clusters: {lambda_len} / { target_clusters }, Time: {(time() - start_22):0.02f}s")
        log_22()
        while lambda_len >= target_clusters:
            i += 1
            clusterMatrix = ClusterSimilarityMatrix.BuildFrom(available_clusters, clusterMatrix)
            alpha1 = clusterMatrix.max_simularity_np()
            log_22()
            if alpha1 < self.aplha1_min:
                break
            else:
                new_available_clusters = MergeClusters(alpha1, clusterMatrix, available_clusters)
        
            if len(new_available_clusters) < target_clusters:
                break
            else:
                available_clusters = new_available_clusters
                lambda_len = len(available_clusters)
        
        self.log("Building simularity matrix (step 2.3)...")
        memberMatrix = MemberMatrix.build(available_clusters, all_items)
        similarity_matrix = memberMatrix.BuildSimularityMatrix(available_clusters)
        certain_clusters = [cluster for cluster in tqdm(available_clusters) if similarity_matrix.Cluster_Max(cluster) > alpha2]
        
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
            alpha2 = similarity_matrix.Cluster_Max(candidate_clusters[len(candidate_clusters)-1])
        
        #calculate new membership and simularity based on candidate and non-candidate
        candidate_memberMatrix = MemberMatrix.build(candidate_clusters, all_items)
        non_candidate_memberMatrix = MemberMatrix.build(non_candidate_clusters, all_items)
        similarity_matrix = MemberSimularityMatrix.RefineTo(similarity_matrix, candidate_clusters)
        
        partition = self.assign_item_to_one_cluster(gamma, alpha2,\
            candidate_clusters, similarity_matrix, candidate_memberMatrix, non_candidate_memberMatrix)
        
        self.log(f"Finished in time {(time() - start_time):0.02f}s")
        
        return partition
    
    
    def assign_item_to_one_cluster(self, gamma: PartitionSet, alpha2: float, candidate_clusters: List[Cluster],\
        similarity_matrix: MemberSimularityMatrix, membership_matrix: MemberMatrix, non_candidate_membership_matrix) -> Partition:
        
        all_items = gamma.get_all_elements()
        
        self.log("Classifying object certainty...")
        totally_certain_lst, certain_lst, uncertain_lst, totally_uncertain_lst \
            = self.identify_object_certainy(all_items, similarity_matrix, alpha2)
        
        self.log("\nFound: ")
        self.log(f"Totally certain objects: {len(totally_certain_lst)}")
        self.log(f"Certain objects: {len(certain_lst)}")
        self.log(f"Uncertain objects: {len(uncertain_lst)}")
        self.log(f"Totally uncertain objects: {len(totally_uncertain_lst)}")
        self.log('\n')
        
        
        totally_uncertain_map, lost_items = {}, []
        if len(totally_uncertain_lst) > 0:
            self.log("Reclassifying totally uncertain items...")
            totally_certain_lst_2, certain_lst_2, uncertain_lst_2, totally_uncertain_map, lost_items =\
                self.reidentify_totally_uncertain_item(totally_uncertain_lst, candidate_clusters,\
                non_candidate_membership_matrix, similarity_matrix, gamma, alpha2)
            
            
            totally_certain_lst += totally_certain_lst_2
            certain_lst += certain_lst_2
            uncertain_lst += uncertain_lst_2
            
            self.log("\nObject certainty after totally uncertain items have been reclassified: ")
            self.log(f"Totally certain objects: {len(totally_certain_lst)}")
            self.log(f"Certain objects: {len(certain_lst)}")
            self.log(f"Uncertain objects: {len(uncertain_lst)}")
            self.log(f"Remaining totally uncertain objects: {len(totally_uncertain_map) + len(lost_items)}")
            self.log('\n')
            
        else:
            self.log("no totally uncertain objects found, skipping reidentification step")
        
        
        self.log("Assigning totally certain objects...")
        candidate_clusters = self.assign_certains_objects(totally_certain_lst, candidate_clusters)
        
        self.log("Assign certain objects...")
        candidate_clusters = self.assign_certains_objects(certain_lst, candidate_clusters)
        
        self.log("Assign uncertain objects")
        candidate_clusters = self.assign_uncertain_objects(uncertain_lst, candidate_clusters, similarity_matrix, membership_matrix)
    
        #Handling lost items. 
        candidate_clusters = self.assign_lost_objects(candidate_clusters, totally_uncertain_map, lost_items)    
                        
        return self.build_final_partition(gamma, candidate_clusters)
    
    
    #Returns tuple of (Totally certain items, certain items, uncertain items Association map, lost items )
    def reidentify_totally_uncertain_item(self, totally_uncertain_item_lst: List, candidate_clusters: List[Cluster],\
        non_can_membership: MemberMatrix, simularity_matrix: MemberSimularityMatrix, gamma: PartitionSet, alpha2: float) -> \
            Tuple[List, List, List, Dict[object, object], List]:
        
        
        self.log("Building coassociation matrix...")
        coassosiation_matrix = CoAssosiationMatrix.build(gamma)
        
        def recalculate_simularity(item: object) -> None:
            Assert.assert_item_in_list(totally_uncertain_item_lst, item)
            for cluster in candidate_clusters:
                v = non_can_membership.average_common_neighbors(coassosiation_matrix, item, cluster)
                simularity_matrix[item, cluster] = v
            return None
        
        self.log('Recalculating simularity using common neighbors...')
        for totally_uncertain_item in tqdm(totally_uncertain_item_lst):
            recalculate_simularity(totally_uncertain_item)
        
        self.log("Rerunning clasification on totally uncertain objects...")
        tot_certain, certain, uncertain, tot_uncertain = self.identify_object_certainy(\
            totally_uncertain_item_lst, simularity_matrix, alpha2)
        tot_uncertain_map, lost_item_lst = {}, []
        
        if len(tot_uncertain) != 0:
            self.log(f'{len(tot_uncertain)} items could not be reidentified, trying to place with most associated item...')
            for uncertain_obj in tot_uncertain:
                #this can contain circular references, watch out
                closest_associate = coassosiation_matrix.FindAssociatedItem(uncertain_obj)
                if closest_associate is not None:
                    tot_uncertain_map[uncertain_obj] = closest_associate
                else:
                    lost_item_lst.append(uncertain_obj)
            
            if len(lost_item_lst) != 0:
                self.log(f'{len(tot_uncertain)} have no associations, placing in isolated clusters')
                #this is done later, as the last thing being calculated
                
        
        return tot_certain, certain, uncertain, tot_uncertain_map, lost_item_lst
        
    
    def assign_lost_objects(self, candidate_clusters: List[Cluster], assosiation_map: Dict[object, object],\
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
    
        
            
    def identify_object_certainy(self, items: List, similarity_matrix: MemberSimularityMatrix, alpha2: float) \
        -> Tuple[List[Tuple[object, Cluster]], List[Tuple[object, Cluster]], List[object], List[object]]:
        totally_certain_lst, certain_lst, uncertain_lst, totally_uncertain_lst  = [], [], [], []
    
        # np.savetxt("matrix.txt", similarity_matrix.matrix, fmt='%f', newline='\n\n')
        # np.savetxt("matrix_sum_0.txt", similarity_matrix.matrix.sum(axis=0), delimiter=',', fmt='%f', newline='\n\n')
        # np.savetxt("matrix_sum_1.txt", similarity_matrix.matrix.sum(axis=1), delimiter=',', fmt='%f', newline='\n\n')
        # np.savetxt("matrix_max_0.txt", similarity_matrix.matrix.max(axis=0), delimiter=',', fmt='%f', newline='\n\n')
        # np.savetxt("matrix_max_1.txt", similarity_matrix.matrix.max(axis=1), delimiter=',', fmt='%f', newline='\n\n')
        reverse_cluster_index_map = {value: key for key, value in similarity_matrix.cluster_index_map.items()}
        reverse_item_index_map = {value: key for key, value in similarity_matrix.item_index_map.items()}

        parameters = [(similarity_matrix.item_index_map[item], similarity_matrix.matrix) for item in items ]
        
        result = None
        with Pool(self.thread_count) as p:
            result = list(tqdm(p.imap(__partial_cluster_certainty_degree__, parameters), total=len(parameters)))
            
        for item_id, similarity, cluster_index in result:
            item_cluster_tuple = (reverse_item_index_map[item_id], reverse_cluster_index_map[cluster_index])
            
            if similarity >= 1:
                totally_certain_lst.append(item_cluster_tuple)
            elif similarity > alpha2 and similarity < 1:
                certain_lst.append(item_cluster_tuple)
            elif similarity <= alpha2 and similarity > 0:
                uncertain_lst.append(item_cluster_tuple[0]) #only add item
            elif similarity == 0:
                totally_uncertain_lst.append(item_cluster_tuple[0]) #only add item
            else:
                raise Exception("something went wrong, simularity value outside bound [0, 1]")
        return totally_certain_lst, certain_lst, uncertain_lst, totally_uncertain_lst
        
    def build_final_partition(self, gamma: PartitionSet, candidate_clusters: List[Cluster]):
        partition = Partition(list(gamma.get_all_elements().keys()))
        
        self.log("Building final partition from candidate clusters...")
        for cluster_index in tqdm(range(len(candidate_clusters))):
            cluster = candidate_clusters[cluster_index]
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



#returns (item_id, s_x, cluster_index)
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