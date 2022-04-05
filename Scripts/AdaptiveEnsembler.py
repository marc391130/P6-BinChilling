from multiprocessing import cpu_count, Pool
from typing import Callable, Dict, List, Tuple
from unittest import result
from Cluster import Cluster, Contig, Partition, PartitionSet
import numpy as np
from tqdm import tqdm
from itertools import islice
import Assertions as Assert
from BitMatrix import BitMatrix
from sys import maxsize as MAXSIZE
from time import time
from MemberSimularityMatrix import MemberSimularityMatrix
from Logger import Logger

from ClusterSimilarityMatrix import ClusterSimilarityMatrix

THREAD_COUNT = min(cpu_count(), 8) 

MERGED_CLUSTER = -1
TOTALLY_CERTAIN, CERTAIN, UNCERTAIN = 3, 2, 1

class Ensembler:
    def ensemble(self, data) -> None:
        pass

class AdaptiveClusterEnsembler(Ensembler):
    def __init__(self, initial_alpha1_thredshold: float = 0.8, initial_delta_aplha: float = 0.1,\
        alpha1_min: float = 0.6, alpha2: float = 0.2, logger: Logger = None):
        self.alpha1_thredshold = initial_alpha1_thredshold
        self.delta_alpha = initial_delta_aplha
        self.aplha1_min = alpha1_min
        self.alpha2 = alpha2
        self.bit_matrix = None
        self.logger = logger                    

    def __stage2__(self, gamma: PartitionSet) -> Tuple[List[Cluster], List[Cluster], Dict[Cluster, int]]:
        initial_clusters = gamma.ClusterToPartitionMap()
        target_clusters = int(gamma.mean_cluster())
        merged_clusters = None
        lambda_len = 0


        print("Started 2.1")
        while True:
            merged_clusters = self.__merge_cls__(self.bit_matrix.similarity_matrix, initial_clusters)
            lambda_len = len(merged_clusters)
            print(f"lambda: {lambda_len}")
            if lambda_len >= target_clusters:
                break
            else:
                self.alpha1_thredshold += self.delta_alpha

        print("Started 2.2")
        i, start_22 = 0, time()
        def log_22() -> None:
            print(f"iteration {i:05d}, alpha1: {self.alpha1_thredshold:.4f} / { self.aplha1_min }, clusters: {lambda_len} / { target_clusters }, Time: {(time() - start_22):0.02f}s")
        log_22()
        while lambda_len >= target_clusters:
            i += 1

            self.alpha1_thredshold = self.bit_matrix.similarity_matrix.max_similarity(list(merged_clusters.keys()))
            log_22()
            if self.alpha1_thredshold < self.aplha1_min:
                break
            else:
                merged_clusters = self.__merge_cls__(self.bit_matrix.similarity_matrix, merged_clusters, True)
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

        #reinit membership matrix
        self.bit_matrix.membership_matrix.initialize(list(cluster_dct.keys()))

        totally_certain_objects_lst, certain_objects_lst, uncertain_objects_lst \
            = self.__find_objects_for_totally_certain_certain_and_uncertain_objects__(list(all_elements.keys()), cluster_dct)
        
        for item in totally_certain_objects_lst:
            Assert.assert_item_in_list(all_elements, item)
            # contains = False
            # for cluster in filtered_cluster_dct.keys():
            #     if item in cluster:
            #         contains = True
            # if not contains:
            #     raise Exception("UGGA BUGGA")
            l = list(filter(lambda cluster: item in cluster, filtered_cluster_dct.keys()))
            if not any(l):
                print("uncertain objects")
                for uncertain_object in totally_uncertain_objects_lst:
                    print(f"item {uncertain_object} with name {uncertain_object.name} is totally uncertain")
                raise Exception(f"item {item} with name {item.name} is ded, len of filter: {len(l)}")
    


        print("Totally certain objects:", len(totally_certain_objects_lst), "Certain Objects:", \
         len(certain_objects_lst), "Uncertain objects:", len(uncertain_objects_lst), "Total elements:", len(all_elements),\
              "Calc total:", str(len(totally_certain_objects_lst) + len(certain_objects_lst) + len(uncertain_objects_lst)))

        # Calculate the quality of candidate clusters

        candidate_clusters_of_totally_certain_objects_dct = {item: [c for c in candidate_clusters if item in c] \
            for item in totally_certain_objects_lst}

        candidate_clusters_of_certain_objects_dct = {item: [c for c in candidate_clusters if item in c] \
            for item in tqdm(certain_objects_lst)}

        print("Started 3.4.1: Finding all uncertain objects in candidate clusters ...")
        # uncertain_objects_lst = self.__find_objects_membership_predicate_for_all__(all_elements, lambda x: x < self.alpha2, filtered_cluster_dct)
        print("Uncertain Objects:", len(uncertain_objects_lst))

        # candidate_clusters_of_uncertain_objects_dct = {item: [c for c in candidate_clusters if item in c] \
        #     for item in uncertain_objects_lst}
        
        # quality_cluster_set = set()
        # for cluster_lst in candidate_clusters_of_uncertain_objects_dct.values():
        #     for cluster in cluster_lst:
        #         if cluster not in quality_cluster_set:
        #             quality_cluster_set.add(cluster)
        
        # quality_dct = {cluster: self.__quality_measure__(cluster, cluster_dct) for cluster in tqdm(quality_cluster_set)}


        print("Evaluating totally certain objects ...")
        self.__assign_certain_objects__(candidate_clusters_of_totally_certain_objects_dct, filtered_cluster_dct)
        print("Evaluating Certain objets ...")
        self.__assign_certain_objects__(candidate_clusters_of_certain_objects_dct, filtered_cluster_dct)
        
        print("Started 3.4.2: Calculating Quality of candidate clusters ...")
        if len(uncertain_objects_lst) > 0:
            quality_dct = {cluster: self.__quality_measure__(cluster, cluster_dct) for cluster in tqdm(candidate_clusters) if len(cluster) != 0}
            
            print("Evaluating uncertain objects...")
            self.__recalculate_and_append_objects__(quality_dct, uncertain_objects_lst, filtered_cluster_dct)
        else:
            print("no uncertain objects to evalutate, skipping step...")
        Assert.assert_only_one_of_each_element_in_cluster_list(all_elements, candidate_clusters)
        #ASSERT.assert_all_elements_are_in_cluster_lst(all_elements, candidate_clusters)

        print("Done")

        return candidate_clusters

    def __assign_certain_objects__(self, certain_objects_map: Dict[object, List[Cluster]], cluster_dct: Dict[Cluster, int]) -> None:
        
        for item, cluster_lst in tqdm(certain_objects_map.items()):
            max_value, max_cluster = 0, None
            Assert.assert_list_nonempty(cluster_lst)
            for cluster in cluster_lst: 
                Assert.assert_list_nonempty(cluster)
                similarity = self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)
                
                if similarity >= max_value:
                    max_value = similarity
                    max_cluster = cluster
                    if similarity >= 1:
                        break
            #end of loop
            if max_cluster is None:
                Assert.assert_not_none(max_cluster)
            for cluster in cluster_lst:
                cluster.remove(item)
                
            max_cluster.append(item)

        
    def __recalculate_and_append_objects__(self, quality_dct: Dict[Cluster, float], uncertain_objects: List, cluster_dct: Dict[Cluster, int]) -> None:
        
        def __partial_quality_measure__(quality_dict: Dict[Cluster, float], cluster: Cluster, item: any, cluster_dct: Dict[Cluster, int]):
            #oldValue is average itemValue / len(cluster). Now that item is added, calculate the contribution of the item.
            #these two are combined by finding the total contribution of oldValue, adding the itemValue and dividing with the length of the cluster
            #Value of oldValue = [total_oldValue / len(cluster - 1)], so we can get total contribution of oldValue by oldValue by:
            #[total_oldValue / len(cluster-1)] * len(cluster-1)
            Assert.assert_item_not_in_collection(cluster, item)
            Assert.assert_list_nonempty(cluster)
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

        for item in tqdm(uncertain_objects):
            min_change, min_change_cluster, new_cluster_quality = MAXSIZE, None, 0

            for cluster in cluster_dct.keys():
                contains_item = item in cluster
                change = 0 if contains_item else __partial_quality_measure__(\
                    quality_dct, cluster, item, cluster_dct)
                new_quality = quality_dct[cluster] + change

                if abs(change) < min_change or change == 0:
                    min_change = change
                    min_change_cluster = cluster
                    new_cluster_quality = new_quality
                    if change == 0:
                        break
            
            if min_change_cluster is not None:
                for cluster in cluster_dct.keys():
                    cluster.remove(item)
                min_change_cluster.append(item)
                quality_dct[min_change_cluster] = new_cluster_quality
        return None

    def __find_objects_for_totally_certain_certain_and_uncertain_objects__(self, elements: List, cluster_dct: Dict[Cluster, int]) -> Tuple[List, List, List]:
        totally_certain_lst, certain_lst, uncertain_lst  = [], [], []
        membership_simularity_matrix = self.bit_matrix.membership_matrix\
            .BuildSimularityMatrix(list(cluster_dct.keys()))
    
        np.savetxt("matrix.txt", membership_simularity_matrix.matrix, fmt='%f', newline='\n\n')
        np.savetxt("matrix_0.txt", membership_simularity_matrix.matrix.sum(axis=0), delimiter=',', fmt='%f', newline='\n\n')
        np.savetxt("matrix_1.txt", membership_simularity_matrix.matrix.sum(axis=1), delimiter=',', fmt='%f', newline='\n\n')
        np.savetxt("matrix_max_0.txt", membership_simularity_matrix.matrix.max(axis=0), delimiter=',', fmt='%f', newline='\n\n')
        np.savetxt("matrix_max_1.txt", membership_simularity_matrix.matrix.max(axis=1), delimiter=',', fmt='%f', newline='\n\n')

        parameters, item_map = [], {} #id -> item
        for i in range(len(elements)):
            item = elements[i]
            tup = (item, i, cluster_dct, membership_simularity_matrix, self.alpha2  )
            item_map[i] = item
            parameters.append( tup )
        
        with Pool(THREAD_COUNT) as p:
            r = list(tqdm(p.imap(__partial_cluster_certainty_degree__, parameters), total=len(parameters)))
            for id, certainty in r:
                item = item_map[id]
                if certainty == TOTALLY_CERTAIN:
                    totally_certain_lst.append(item)
                elif certainty == CERTAIN:
                    certain_lst.append(item)
                elif certainty == UNCERTAIN:
                    uncertain_lst.append(item)
                else:
                    raise Exception("something went wrong")
        
        # for item in tqdm(elements):
        #     uncertain = True
        #     for cluster in cluster_dct.keys():

        #         similarity = self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)

        #         if similarity >= 1:
        #             totally_certain_lst.append(item)
        #             uncertain = False
        #             break
        #         elif similarity > self.alpha2 and similarity < 1:
        #             certain_lst.append(item)
        #             uncertain = False
        #             break

        #     if uncertain:
        #         uncertain_lst.append(item)
                    
        return totally_certain_lst, certain_lst, uncertain_lst
    
    def __find_objects_membership_predicate_for_all__(self, elements: List, predicate : Callable[[float], bool], cluster_dct: Dict[Cluster, int]) -> List:
        objects_lst = []
        for item in tqdm(elements):
            uncertain = True
            for cluster, _ in cluster_dct.items():
                if item in cluster:
                    uncertain = False
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
        for uncertain_object, cluster_lst in tqdm(non_candidate_clusters_of_uncertain_objects_dct.items()):
            max_similarity, max_similarity_candidate_cluster = -1, None
            Assert.assert_list_nonempty(cluster_lst)
            for cluster in cluster_lst:
                for candidate_cluster in candidate_clusters:
                    similarity = self.bit_matrix.similarity_matrix.getEntry(cluster, candidate_cluster)
                    # similarity = gamma.__similarity_measure_cluster__(cluster, candidate_cluster)
                    if similarity >= max_similarity:
                        max_similarity = similarity
                        max_similarity_candidate_cluster = candidate_cluster
            
            Assert.assert_not_none(max_similarity_candidate_cluster)
            if uncertain_object in max_similarity_candidate_cluster:
                raise Exception(f"Item {uncertain_object} is already in cluster!")
            max_similarity_candidate_cluster.append(uncertain_object)

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



    def ensemble(self, gamma: PartitionSet[Contig]) -> None:
        self.bit_matrix = BitMatrix(gamma)
        # print(self.bit_matrix, self.bit_matrix.matrix.shape)
        candidate_clusters, non_candidate_clusters, merged_clusters = self.__stage2__(gamma)
        print(len(candidate_clusters), len(non_candidate_clusters), self.alpha2, gamma.maximal_partition_clusters())
        return self.__stage3__(candidate_clusters, non_candidate_clusters, gamma, merged_clusters)

    def print_to_file(self, file_path: str, clusters: List[Cluster]) -> None:
        with open(file_path, 'w') as file:
            for cluster_idx in range(len(clusters)):
                for item in clusters[cluster_idx]:
                    file.write(f"{cluster_idx}\t{item.name}\n")

    def __is_same_partition__(self, partition_idx1: int, partition_idx2: int) -> bool:
            return partition_idx1 == partition_idx2 and partition_idx1 != MERGED_CLUSTER

    def __merge_cls__(self, similarity_matrix: ClusterSimilarityMatrix, dct_info: Dict[Cluster, int], debug:bool = False) -> Dict[Cluster, int]:
        result_dct = {}
        done_set = set()
        list_info = list(dct_info.items())
        def getCluster(index: int) -> Cluster:
            return list_info[index][0]
        def is_same_partition(i, j) -> bool:
            return self.__is_same_partition__(list_info[i][1], list_info[j][1])

        for i in range(len(list_info)):
            cluster1 = getCluster(i)
            add_to_results = True
            if cluster1 in done_set:
                continue

            for j in range(i, len(list_info)):
                cluster2 = getCluster(j)

                if cluster1 is cluster2 or is_same_partition(i, j) or cluster2 in done_set:
                    continue
                    
                if similarity_matrix[(cluster1, cluster2)] >= self.alpha1_thredshold:
                    new_cluster = Cluster.merge(cluster1, cluster2)
                    result_dct[new_cluster] = MERGED_CLUSTER
                    done_set.add(cluster1)
                    done_set.add(cluster2)
                    self.bit_matrix.add_cluster_to_bit_matrix(new_cluster)
                    add_to_results = False
                    break
            
            if add_to_results:
                result_dct[cluster1] = dct_info[cluster1]
                
        return result_dct


    # def identify_data_certainty(self, gamma : PartitionSet, cluster_dct: Dict[Cluster, int]) -> tuple[List, List, List, List]:
    #     certain_lst, uncertain_lst, totally_certain_lst, totally_uncertain_lst = [], [], [], []
    #     all_items = gamma.get_all_elements()
        
    #     print(f"splitting ({len(all_items)}) items into certainty groups using '{len(cluster_dct)}' clusters...")
    #     for item in tqdm(all_items):  
    #         is_certain, is_uncertain, is_totally_uncertain = False, False, True                       
    #         for cluster in cluster_dct.keys():

    #             similarity = self.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)

    #             if similarity >= 1:
    #                 totally_certain_lst.append(item)
    #                 is_certain, is_uncertain, is_totally_uncertain = False, False, False
    #                 break
    #             elif similarity > self.alpha1_thredshold:
    #                 is_certain, is_uncertain, is_totally_uncertain = True, False, False
    #             elif (not is_certain) and similarity <= self.alpha1_thredshold and similarity != 0:
    #                 is_certain, is_uncertain, is_totally_uncertain = False, True, False
            
    #         if is_certain:
    #             certain_lst.append(item)
    #             continue
    #         if is_uncertain:
    #             uncertain_lst.append(item)
    #             continue
    #         if is_totally_uncertain:
    #             totally_uncertain_lst.append(item)
    #             continue
    #     return totally_certain_lst, certain_lst, uncertain_lst, totally_uncertain_lst
    
    
    
def __partial_cluster_certainty_degree__(\
    tuple : Tuple[object, int,  Dict[Cluster, int], MemberSimularityMatrix, float] ) -> Tuple[object, int]:
    # item, id, cluster_dct, ensembler = tuple
    item, id, cluster_dct, matrix, alpha2 = tuple
    
    certainty_degree = UNCERTAIN
    for cluster in cluster_dct.keys():

        similarity = matrix.GetEntry(cluster, item)
        # similarity = ensembler.bit_matrix.membership_similarity_measure(item, cluster, cluster_dct)

        if similarity >= 1:
            certainty_degree = TOTALLY_CERTAIN
            break
        elif similarity > alpha2 and similarity < 1:
            certainty_degree = CERTAIN
            
            
    return (id, certainty_degree)