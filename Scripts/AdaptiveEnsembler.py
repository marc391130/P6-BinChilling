from ast import Assert
from typing import Dict, List, Tuple
from unittest import result
from Cluster import Cluster, Contig, Partition, PartitionSet
import numpy as np
from tqdm import tqdm
from itertools import islice
import Assertions as ASSERT
from BitMatrix import BitMatrix
from sys import maxsize as MAXSIZE

MERGED_CLUSTER = -1

class Ensembler:
    def ensemble(self, data) -> None:
        pass

class AdaptiveClusterEnsembler(Ensembler):
    def __init__(self, initial_alpha1_thredshold: float = 0.8, initial_delta_aplha: float = 0.1, alpha1_min: float = 0.6, alpha2: float = 0.2):
        self.alpha1_thredshold = initial_alpha1_thredshold
        self.delta_alpha = initial_delta_aplha
        self.aplha1_min = alpha1_min
        self.alpha2 = alpha2
        self.bit_matrix = None


    def generate_consensus_clusters(self, gamma: PartitionSet, bit_matrix: np.matrix) -> List[Cluster]:
        dct_info = self.__map_clusters__(gamma) # dict [Cluster, partiton_idx]

        target_clusters = gamma.maximal_partition_clusters()
        while True:
            result_dct = self.__merge_similar_clusters__(gamma, dct_info)
            lambda_clusters = len(result_dct)
            if lambda_clusters > target_clusters:
                # Lower alpha1
                pass
            elif lambda_clusters < target_clusters:
                pass
            else:
                return result_dct                
            
    def __map_clusters__(self, gamma: PartitionSet) -> Dict[Cluster, int]:
        dct_info = {}
        for partition_idx in range(len(gamma)):
            for key, value in gamma[partition_idx].items():
                dct_info[value] = partition_idx
        return dct_info

    def __stage2__(self, gamma: PartitionSet) -> Tuple[List[Cluster], List[Cluster], Dict[Cluster, int]]:
        initial_clusters = self.__map_clusters__(gamma)
        target_clusters = int(gamma.mean_cluster())
        similarity_matrix = self.__similarity_matrix__(gamma, initial_clusters)
        merged_clusters = None
        lambda_len = 0


        print("Started 2.1")
        while True:
            merged_clusters = self.__merge_cls__(similarity_matrix, initial_clusters)
            lambda_len = len(merged_clusters)
            if lambda_len >= target_clusters:
                break
            else:
                self.alpha1_thredshold += self.delta_alpha

        print("Started 2.2")
        i = 0
        while lambda_len >= target_clusters:

            if i % 1 == 0:
                print("Current Loop", i, self.alpha1_thredshold, lambda_len, target_clusters)
            i += 1

            similarity_matrix = self.__similarity_matrix__(gamma, merged_clusters)
            self.alpha1_thredshold = similarity_matrix.max()
            if self.alpha1_thredshold < self.aplha1_min:
                break
            else:
                merged_clusters = self.__merge_cls__(similarity_matrix, merged_clusters, True)
                lambda_len = len(merged_clusters)
            


        print("Started 2.3")
        certain_clusters = self.__find_all_clusters_with_atleast_one_certain_cluster__(merged_clusters)

        candidate_clusters = []
        non_candidate_clusters = []

        if len(certain_clusters) == target_clusters:
            candidate_clusters = certain_clusters
            non_candidate_clusters = [k for k,v in merged_clusters.items() if k not in certain_clusters]
        else:
            cluster_certainty_lst = [(k, v, self.__mean_membership_similarity_measure__(k, merged_clusters)) for k, v in merged_clusters.items()]
            cluster_certainty_lst = list(sorted(cluster_certainty_lst, key = lambda item: item[2], reverse=True))
            candidate_clusters = [x[0] for x in islice(cluster_certainty_lst, target_clusters)]
            self.alpha2 = self.__max_membership_similarity_measure__(candidate_clusters[len(candidate_clusters) - 1], merged_clusters)
            non_candidate_clusters = [x[0] for x in islice(cluster_certainty_lst, target_clusters, None)]

        return (candidate_clusters, non_candidate_clusters, merged_clusters)

    def __stage3__(self, candidate_clusters: List[Cluster], non_candidate_clusters: List[Cluster], gamma: PartitionSet, cluster_dct: Dict[Cluster, int]) -> List[Cluster]:
        # Find all totally uncertain objects in candicate clusters 
        # (Aka The objects that are not in any candidate cluster)
        
        print("Started 3.1: Finding totally uncertain objects ...")
        filtered_cluster_dct = dict(filter(lambda elem: elem[0] in candidate_clusters, cluster_dct.items()))
        all_elements = gamma.get_all_elements()

        totally_uncertain_objects_lst = self.__find_objects_membership_predicate_for_all__(all_elements, lambda x: x <= 0, filtered_cluster_dct)

        print("Totally uncertain objects:", len(totally_uncertain_objects_lst))
        non_candidate_clusters_of_uncertain_objects_dct = {item: [c for c in non_candidate_clusters if item in c] \
            for item in totally_uncertain_objects_lst}
        
        # Add uncertain objects to the candidate cluster that is most similair to the non-candidate cluster with the uncertain object'
        print("Started 3.2: Adding totally uncertain objects to candidate clusters from non-candidate clusters ...")
        self.__add_totally_uncertain_objects_to_candidate_clusters__(non_candidate_clusters_of_uncertain_objects_dct, candidate_clusters, gamma)

        # Find totally certain objects and certain objects in candidate clusters
        print("Started 3.3.1: Finding all totally certain and certain objects ...")

        totally_certain_objects_lst, certain_objects_lst, uncertain_objects_lst = self.__find_objects_for_totally_certain_certain_and_uncertain_objects__(all_elements, filtered_cluster_dct)

        # totally_certain_objects_lst = self.__find_objects_membership_predicate__(all_elements, lambda x: (x >= 1), filtered_cluster_dct)
        # certain_objects_lst = self.__find_objects_membership_predicate__(all_elements, lambda x: (x >= self.alpha2 and x < 1), filtered_cluster_dct)

        print("Totally certain objects:", len(totally_certain_objects_lst), "Certain Objects:", \
         len(certain_objects_lst), "Uncertain objects:", len(uncertain_objects_lst), "Total elements:", len(all_elements),\
              "Calc total:", str(len(totally_certain_objects_lst) + len(certain_objects_lst) + len(uncertain_objects_lst)))

        # Calculate the quality of candidate clusters
        print("Started 3.3.2: Calculating Quality of candidate clusters ...")
        
        quality_dct = {cluster: self.__quality_measure__(cluster, cluster_dct) for cluster in tqdm(candidate_clusters)}

        candidate_clusters_of_totally_certain_objects_dct = {item: [c for c in candidate_clusters if item in c] \
            for item in totally_certain_objects_lst}

        candidate_clusters_of_certain_objects_dct = {item: [c for c in candidate_clusters if item in c] \
            for item in tqdm(certain_objects_lst)}

        print("Started 3.4: Finding all uncertain objects in candidate clusters ...")
        # uncertain_objects_lst = self.__find_objects_membership_predicate_for_all__(all_elements, lambda x: x < self.alpha2, filtered_cluster_dct)
        print("Uncertain Objects:", len(uncertain_objects_lst))

        candidate_clusters_of_uncertain_objects_dct = {item: [c for c in candidate_clusters if item in c] \
            for item in uncertain_objects_lst}


        print("Evaluating totally certain objects ...")
        self.__recalculate_and_append_objects__(quality_dct, candidate_clusters_of_totally_certain_objects_dct, filtered_cluster_dct)
        print("Evaluating Certain objets ...")
        self.__recalculate_and_append_objects__(quality_dct, candidate_clusters_of_certain_objects_dct, filtered_cluster_dct)
        print("Evaluating uncertain objects ...")
        self.__recalculate_and_append_objects__(quality_dct, candidate_clusters_of_uncertain_objects_dct, filtered_cluster_dct)
        ASSERT.assert_only_one_of_each_element_in_cluster_list(all_elements, candidate_clusters)
        #ASSERT.assert_all_elements_are_in_cluster_lst(all_elements, candidate_clusters)

        print("Done")

        return candidate_clusters

        
    def __recalculate_and_append_objects__(self, quality_dct: Dict[Cluster, float], obj_list_cluster_dct: Dict[object, List[Cluster]], cluster_dct: Dict[Cluster, int]) -> None:
        
        def __partial_quality_measure__(self, quality_dict: Dict[Cluster, float], cluster: Cluster, item: any, cluster_dct: Dict[Cluster, int]):
            #oldValue is average itemValue / len(cluster). Now that item is added, calculate the contribution of the item.
            #these two are combined by finding the total contribution of oldValue, adding the itemValue and dividing with the length of the cluster
            #Value of oldValue = [total_oldValue / len(cluster - 1)], so we can get total contribution of oldValue by oldValue by:
            #[total_oldValue / len(cluster-1)] * len(cluster-1)
            ASSERT.assert_min_list_len(cluster, 1)
            oldValue = quality_dict[cluster]
            total_oldValue = oldValue * (len(cluster))
            itemValue = 0
            try:
                cluster.append(item)
                itemValue = pow((self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct) \
                    - self.__mean_membership_similarity_measure__(cluster, cluster_dct)), 2)
            finally:
                cluster.remove(item)
            
            return (total_oldValue + itemValue) / (len(cluster)+1)

        for item, cluster_lst in tqdm(obj_list_cluster_dct.items()):
            min_change, min_change_cluster, new_cluster_quality = MAXSIZE, None, 0

            for cluster in cluster_lst:
                contains_item = item in cluster
                change = 0 if contains_item else __partial_quality_measure__(\
                    self, quality_dct, cluster, item, cluster_dct)
                new_quality = quality_dct[cluster] + change

                if abs(change) < min_change or change == 0:
                    min_change = change
                    min_change_cluster = cluster
                    new_cluster_quality = new_quality
                    if change == 0:
                        break
            
            if min_change_cluster is not None:
                for cluster in cluster_lst:
                    cluster.remove(item)
                min_change_cluster.append(item)
                quality_dct[min_change_cluster] = new_cluster_quality
        return None
    
    def __find_objects_membership_predicate__(self, elements: List, predicate, cluster_dct: Dict[Cluster, int]) -> List:
        objects_lst = []
        check_set = set()

        for item in tqdm(elements):
            for cluster, _ in cluster_dct.items():
                if predicate(self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)):
                    if item not in check_set:
                        objects_lst.append(item)
                        check_set.add(item)
                        break
        return objects_lst

    def __find_objects_for_totally_certain_certain_and_uncertain_objects__(self, elements: List, cluster_dct: Dict[Cluster, int]) -> Tuple[List, List, List]:
        totally_certain_lst = []
        certain_lst = []
        uncertain_lst = []
        
        for item in tqdm(elements):
            uncertain = True
            for cluster in cluster_dct.keys():

                similarity = self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)

                if similarity >= 1:
                    totally_certain_lst.append(item)
                    uncertain = False
                    break
                elif similarity > self.alpha2 and similarity < 1:
                    certain_lst.append(item)
                    uncertain = False
                    break

            if uncertain:
                uncertain_lst.append(item)
                    
        return totally_certain_lst, certain_lst, uncertain_lst
    
    def __find_objects_membership_predicate_for_all__(self, elements: List, predicate, cluster_dct: Dict[Cluster, int]) -> List:
        objects_lst = []
        not_uncertain = 0
        for item in tqdm(elements):
            uncertain = True
            for cluster, _ in cluster_dct.items():
                if not predicate(self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)):
                    uncertain = False
                    not_uncertain += 1
                    break
            if uncertain:
                objects_lst.append(item)
                        
        return objects_lst

    def __quality_measure__(self, cluster: Cluster, cluster_dct: Dict[Cluster, int]) -> float:
        if len(cluster) == 0:
            raise Exception("Trying to get Quality of cluster that is empty ;-(")

        result = 0
        for item in cluster:
            result += pow((self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct) - self.__mean_membership_similarity_measure__(cluster, cluster_dct)), 2)
        return result / len(cluster)


    def __add_totally_uncertain_objects_to_candidate_clusters__(self, non_candidate_clusters_of_uncertain_objects_dct: Dict[str, List[Cluster]], candidate_clusters: List[Cluster], gamma: PartitionSet) -> None:
        failed = 0
        for uncertain_object, cluster_lst in non_candidate_clusters_of_uncertain_objects_dct.items():
            max_similarity = 0
            max_similarity_candidate_cluster = None
            for cluster in cluster_lst:
                for candidate_cluster in candidate_clusters:
                    similarity = gamma.__similarity_measure_cluster__(cluster, candidate_cluster)
                    if similarity > max_similarity:
                        max_similarity = similarity
                        max_similarity_candidate_cluster = candidate_cluster
            
            if max_similarity_candidate_cluster is not None:
                if uncertain_object in max_similarity_candidate_cluster:
                    print(f"Object {str(uncertain_object)} already in cluster!")
                    failed += 1
                else:
                    max_similarity_candidate_cluster.append(uncertain_object)
        print(f"Failed: {failed}")

    def __find_all_clusters_with_atleast_one_certain_cluster__(self, cluster_dct: Dict[Cluster, int]) -> List[Cluster]:
        result = []

        all_merged_cluster = dict(filter(lambda x: x[1] == MERGED_CLUSTER, cluster_dct.items()))
        
        for cluster, partition_idx in all_merged_cluster.items():
            for item in cluster:
                membership_measure = self.bit_matrix.membership_similarity_measure(item, cluster, all_merged_cluster)
                if membership_measure > self.alpha2:
                    result.append(cluster)
                    break

        return result

    def __max_membership_similarity_measure__(self, cluster: Cluster, cluster_dct: Dict[Cluster, int]) -> float:
        result = 0

        for item in cluster:
            result = max(self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct), result)
        
        return result

    def __mean_membership_similarity_measure__(self, cluster: Cluster, cluster_dct: Dict[Cluster, int]) -> float:
        result = 0

        for item in cluster:
            result += self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)
        
        return result / len(cluster)


    def __similarity_matrix__(self, gamma: PartitionSet, cluster_dct: Dict[Cluster, int]) -> np.matrix:
        matrix = np.empty_like(0, shape = (len(cluster_dct), len(cluster_dct))).astype(np.float32)
        matrix.fill(0)
        lst_data = list(cluster_dct.items())
        done_set = set()

        for i in range(len(lst_data)):
            for j in range(len(lst_data)):
                # Also catches if cluster1 == cluster 2, as they are in same partition
                if self.__is_same_partition__(lst_data[i][1], lst_data[j][1]) or i == j:
                    continue

                if lst_data[j][0] in done_set:
                    continue
                
                similarity = gamma.__similarity_measure_cluster__(lst_data[i][0], lst_data[j][0])
                matrix[i,j] = similarity
                matrix[j,i] = similarity

        return matrix

    def __is_same_partition__(self, partition_idx1: int, partition_idx2: int) -> bool:
        return partition_idx1 == partition_idx2 and partition_idx1 != MERGED_CLUSTER
            

    def __merge_similar_clusters__(self, gamma: PartitionSet, dct_info: Dict[Cluster, int], done_dct: set = set()) -> Dict[Cluster, int]:
        similarities_found = 0
        result_dct = {}
        done_dct = set()
        for cluster1, partition_idx1 in dct_info.items():
            if cluster1 in done_dct:
                continue

            for cluster2, partition_idx2 in dct_info.items():
                if cluster1 is cluster2 or (partition_idx1 == partition_idx2 and partition_idx1 != MERGED_CLUSTER) or cluster2 in done_dct:
                    continue
            
                similarity = gamma.__similarity_measure_cluster__(cluster1, cluster2)

                if similarity >= self.alpha1_thredshold:
                    new_cluster = cluster1.merge(cluster2)
                    result_dct[new_cluster] = MERGED_CLUSTER
                    done_dct.add(cluster1)
                    done_dct.add(cluster2)
                    similarities_found += 1
                    break
            
            if cluster1 not in done_dct:
                result_dct[cluster1] = dct_info[cluster1]
        
        if similarities_found == 0:
            return result_dct
        else:
            return self.__merge_similar_clusters__(gamma, result_dct, done_dct)


    def ensemble(self, gamma: PartitionSet[Contig]) -> None:
        self.bit_matrix = BitMatrix(gamma)
        print(self.bit_matrix, self.bit_matrix.matrix.shape)
        candidate_clusters, non_candidate_clusters, merged_clusters = self.__stage2__(gamma)
        print(len(candidate_clusters), len(non_candidate_clusters), self.alpha2, gamma.maximal_partition_clusters())
        return self.__stage3__(candidate_clusters, non_candidate_clusters, gamma, merged_clusters)

    def print_to_file(self, file_path: str, clusters: List[Cluster]) -> None:
        with open(file_path, 'w') as file:
            for cluster_idx in range(len(clusters)):
                for item in clusters[cluster_idx]:
                    file.write(f"{cluster_idx}\t{item.name}\n")

    def __merge_cls__(self, similarity_matrix: np.matrix, dct_info: Dict[Cluster, int], debug:bool = False) -> Dict[Cluster, int]:
        result_dct = {}
        done_set = set()
        list_info = list(dct_info.items())
        merges = 0

        for i in range(len(list_info)):
            add_to_results = True
            if list_info[i][0] in done_set:
                continue

            for j in range(i, len(list_info)):

                if list_info[i][0] is list_info[j][0] or self.__is_same_partition__(list_info[i][1], list_info[j][1]) or list_info[j][0] in done_set:
                    continue
                    
                if similarity_matrix[i][j] >= self.alpha1_thredshold:
                    merges += 1
                    new_cluster = list_info[i][0].merge(list_info[j][0])
                    result_dct[new_cluster] = MERGED_CLUSTER
                    done_set.add(list_info[i][0])
                    done_set.add(list_info[j][0])
                    add_to_results = False
                    break
            
            if add_to_results:
                key = list_info[i][0]
                result_dct[key] = dct_info[key]
                self.bit_matrix.add_cluster_to_bit_matrix(key)

        return result_dct

