import numpy as np
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from typing import Tuple, Dict, List
import Assertions as ASSERT
from ClusterSimilarityMatrix import ClusterSimilarityMatrix
import formulas as formula
import Constants as Constant

class BitMatrix:
    def __init__(self, gamma: PartitionSet):
        self.matrix = None #data
        self.__gamma__ = gamma
        self.element_index_map = None #columns
        self.cluster_index_map = None #rows
        self.similarity_matrix = None
        self.bit_matrix_transform(gamma)


    def bit_matrix_transform(self, gamma: PartitionSet) -> None:
        clusters = gamma.get_all_clusters()
        self.element_index_map = gamma.get_all_elements()
        
        self.cluster_index_map = self.__setup_cluster_dct__(clusters)
        self.matrix = self.__setup_bit_matrix__(self.cluster_index_map, self.element_index_map)
        self.similarity_matrix = self.__build_similarity_matrix__(gamma, gamma.ClusterToPartitionMap())

        self.__update_matrix__()
        
        
    def __build_similarity_matrix__(self, gamma : PartitionSet, cluster_dct: Dict[Cluster, int]) -> ClusterSimilarityMatrix:
        size = len(self.cluster_index_map)
        matrix = np.full(shape = (size, size), fill_value=np.nan)
        total_elements = gamma.__total_elements__()
        
        print("calculating similarity matrix...")
        for cluster1, p1 in tqdm(cluster_dct.items()):
            index_1 = self.cluster_index_map[cluster1]
            for cluster2, p2 in cluster_dct.items():                
                index_2 = self.cluster_index_map[cluster2]
                if cluster1 is cluster2:
                    matrix[index_1, index_2] = np.nan
                    continue
                
                if matrix[index_1, index_2] == np.nan:
                    continue #skip if already calculated
                
                value = formula.cluster_simularity(cluster1, cluster2, len(cluster1.intersection(cluster2)), total_elements) \
                    if p1 != p2 else np.nan
                
                matrix[index_1, index_2] = value
                matrix[index_2, index_1] = value
                
        return ClusterSimilarityMatrix(matrix, self.cluster_index_map)
    
    def __build_similarity_row__(self, cluster: Cluster, partition: int) -> Dict[int, float or np.nan]:
        result = {}
        total_elements = self.__gamma__.__total_elements__()
        ASSERT.assert_key_exists(cluster, self.cluster_index_map)
        exception_lst = set(cluster.__get_leaf_clusters__())
        index_1 = self.cluster_index_map[cluster]
        for cluster2, p2 in self.cluster_index_map.items():                
            index_2 = self.cluster_index_map[cluster2]
            if cluster is cluster2 or index_1 == index_2 or cluster2 in exception_lst:
                result[index_2] = np.nan
                continue           
    
            value = formula.cluster_simularity(cluster, cluster2, len(cluster.intersection(cluster2)), total_elements) \
                if partition != p2 or partition != Constant.MERGED_CLUSTER else np.nan
            
            result[index_2] = value    
        return result    
        

    def add_cluster_to_bit_matrix(self, cluster: Cluster) -> None:
        
        ASSERT.assert_key_not_exists(cluster, self.cluster_index_map)
        ASSERT.assert_not_none(self.matrix)
        ASSERT.assert_not_none(self.similarity_matrix)
        # ASSERT.assert_item_not_in_collection(self.cluster_index_map, cluster)
        new_index = len(self.cluster_index_map)
        self.cluster_index_map[cluster] = new_index
        self.matrix = np.append(self.matrix, np.zeros( shape=(len(self.element_index_map), 1) ), axis=1)
        
        similarity_row = self.__build_similarity_row__(cluster, Constant.MERGED_CLUSTER)
        self.similarity_matrix.__expand__(cluster, similarity_row)
        
        self.__update_column__(cluster)
        
        # self.matrix = self.__setup_bit_matrix__(self.cluster_index_map, self.element_index_map)

        # self.__update_matrix__()

    def __update_matrix__(self) -> None:
        for cluster, _ in self.cluster_index_map.items():
            self.__update_column__(cluster)

    def __update_column__(self, cluster: Cluster):
        cluster_idx = self.cluster_index_map[cluster]
        for item in cluster:
                self.matrix[self.element_index_map[item], cluster_idx] = 1

    def calc_membership_in_cluster(self, item, cluster: Cluster) -> int:
        result = 0
        for leaf in cluster.__get_leaf_clusters__():
            result += self.get_entry(item, leaf)
        return result
        
    def membership_similarity_measure(self, item, cluster: Cluster, cluster_dct: Dict[Cluster, int]) -> float:
        # all_merged_cluster = dict(filter(lambda x: x[1] == MERGED_CLUSTER, cluster_dct.items())) if use_filter else cluster_dct
        # MIGHT HAVE FUCKED UP ^^^ BUT IT works maybe. Maybe filter is required, however it fixes bugs to not use.
        max_value, cluster_value = 0, 0

        for cluster__, partition_idx in cluster_dct.items():
            membership_value = cluster__.calc_membership(item)
            max_value = max(membership_value, max_value)
            if cluster__ is cluster:
                cluster_value = membership_value
            
        if max_value == 0:
            #return 0
            raise Exception(f"The item {str(item)} does not exist in the dataset")

        return cluster_value / max_value

    def get_entry(self, item, cluster: Cluster) -> int:
        return self.matrix[self.element_index_map[item], self.cluster_index_map[cluster]]
    
    #object__: first argument is item, second is the cluster
    def __getitem__(self, object__: Tuple[any, Cluster]) -> int:
        item, cluster = object__
        return self.get_entry(item, cluster)

    def __setup_cluster_dct__(self, clusters: List[Cluster]) -> Dict[Cluster, int]:
        result = {}
        for cluster_idx in range(len(clusters)):
            result[clusters[cluster_idx]] = cluster_idx
        return result


    def __setup_bit_matrix__(self, clusters, elements) -> np.matrix:
        column_nr = len(clusters) + 1
        row_nr = len(elements)

        matrix = np.empty_like(0, shape = (row_nr, column_nr))
        matrix.fill(0)
        return matrix

