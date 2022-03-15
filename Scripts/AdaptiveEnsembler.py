from typing import Dict, List, Tuple
from unittest import result
from Cluster import Cluster, Contig, Partition, PartitionSet
import numpy as np
from tqdm import tqdm
import itertools

MERGED_CLUSTER = -1

class Ensembler:
    def ensemble(self, data) -> None:
        pass

class AdaptiveClusterEnsembler(Ensembler):
    def __init__(self, initial_alpha1_thredshold: float = 0.8, initial_delta_aplha: float = 0.1, alpha1_min: float = 0.6):
        self.alpha1_thredshold = initial_alpha1_thredshold
        self.delta_alpha = initial_delta_aplha
        self.aplha1_min = alpha1_min

    def bit_matrix_transform(self, gamma: PartitionSet) -> np.matrix:
        matrix = self.__setup_bit_matrix(gamma)
        
        values_idx = gamma.get_all_elements()
        clusters = gamma.get_all_clusters()
        for cluster_idx in tqdm(range(len(clusters))):
            for item in clusters[cluster_idx]:
                matrix[values_idx[item], cluster_idx] = True
        
        return matrix

    def __setup_bit_matrix(self, gamma: PartitionSet) -> np.matrix:
        column_nr = len(gamma.get_all_clusters())
        row_nr = len(gamma.get_all_elements())

        matrix = np.empty_like(False, shape = (row_nr, column_nr))
        matrix.fill(False)
        return matrix

    def generate_consensus_clusters(self, gamma: PartitionSet, bit_matrix: np.matrix) -> List[Cluster]:
        dct_info = self.__map_clusters__(gamma) # dict [Cluster, partiton_idx]

        target_clusters = gamma.maximal_partition_clusters()
        while True:
            result_dct = self.__merge_similar_clsuters__(gamma, dct_info)
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

    def __stage2__(self, gamma: PartitionSet, bit_matrix: np.matrix) -> Dict[Cluster, int]:
        initial_clusters = self.__map_clusters__(gamma)
        target_clusters = gamma.maximal_partition_clusters()
        similarity_matrix = self.__similarity_matrix__(gamma, initial_clusters)
        merged_clusters = None
        lambda_len = 0

        while True:
            merged_clusters = self.__merge_cls__(similarity_matrix, initial_clusters)
            lambda_len = len(merged_clusters)
            if lambda_len >= target_clusters:
                break
            else:
                self.alpha1_thredshold += self.delta_alpha

        while lambda_len >= target_clusters:
            similarity_matrix = self.__similarity_matrix__(gamma, merged_clusters)
            self.alpha1_thredshold = similarity_matrix.max()
            if self.alpha1_thredshold < self.aplha1_min:
                break
            else:
                merged_clusters = self.__merge_cls__(similarity_matrix, merged_clusters)
                lambda_len = len(merged_clusters)
            
            
    def __similarity_matrix__(self, gamma: PartitionSet, cluster_dct: Dict[Cluster, int]) -> np.matrix:
        matrix = np.empty_like(0, shape = (len(cluster_dct), len(cluster_dct))).astype(np.float32)
        matrix.fill(0)
        done_set = set()
        i = 0

        for cluster1, partition1_idx in cluster_dct.items():
            j = 0

            for cluster2, partition2_idx in cluster_dct.items():
                # Also catches if cluster1 == cluster 2, as they are in same partition
                if self.__is_same_partition__(partition1_idx, partition2_idx) or i == j:
                    continue

                if cluster2 in done_set:
                    matrix[i,j] = matrix[j,i]
                    continue
                
                matrix[i,j] = gamma.__similarity_measure_cluster__(cluster1, cluster2)

                j += 1
            i += 1

        return matrix

    def __is_same_partition__(self, partition_idx1: int, partition_idx2: int) -> bool:
        return partition_idx1 == partition_idx2 and partition_idx1 != MERGED_CLUSTER
            

    def __merge_similar_clsuters__(self, gamma: PartitionSet, dct_info: Dict[Cluster, int], done_dct: set = set()) -> Dict[Cluster, int]:
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
            return self.__merge_similar_clsuters__(gamma, result_dct, done_dct)

    def __membership_similarity__(self, bit_matrix: np.matrix, item, cluster: Cluster, merged_lst: List[Cluster]) -> None:
        pass


    def ensemble(self, gamma: PartitionSet[Contig]) -> None:
        data = self.bit_matrix_transform(gamma)
        print(data, data.shape)
        merged = self.__stage2__(gamma, data)
        print(merged, len(merged))

    def __merge_cls__(self, similarity_matrix: np.matrix, dct_info: Dict[Cluster, int]) -> Dict[Cluster, int]:
        result_dct = {}
        done_set = set()
        list_info = list(dct_info.items())

        for i in range(len(list_info)):
            if list_info[i][0] in done_set:
                continue

            for j in range(len(list_info)):
                if list_info[i][0] is list_info[j][0] or self.__is_same_partition__(list_info[i][1], list_info[j][1]) or list_info[j][0] in done_set:
                    continue

                if similarity_matrix[i,j] >= self.alpha1_thredshold:
                    new_cluster = list_info[i][0].merge(list_info[j][0])
                    result_dct[new_cluster] = MERGED_CLUSTER
                    done_set.add(list_info[i][0])
                    done_set.add(list_info[j][0])
                    break
            
            if list_info[i][0] not in done_set:
                result_dct[list_info[i][0]] = dct_info[list_info[i][0]]

        return result_dct

