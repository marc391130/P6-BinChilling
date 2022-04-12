from __future__ import annotations
import numpy as np
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from typing import Tuple, Dict, List, Callable
import Assertions as Assert 
from math import sqrt

class ClusterSimilarityMatrix():
    def __init__(self, matrix: np.matrix, index_map: Dict[Cluster, int], item_count: int) -> None:
        self.matrix = matrix
        self.index_map = index_map
        self.__item_count__ = item_count
        
    @staticmethod
    def Build(clusters: List[Cluster], item_count: int) -> ClusterSimilarityMatrix:
        cluster_index_map = {clusters[index]: index for index in range(len(clusters)) }
        
        size = len(cluster_index_map)
        matrix = np.full(shape = (size, size), fill_value=np.NINF, dtype=np.float32)
        
        #         matrix[matrix_index_2, matrix_index_1] = value
        def calculate(c1: Cluster, c2: Cluster) -> float:
            return ClusterSimilarityMatrix.cluster_simularity(c1, c2, len(c1.intersection(c2)), item_count)
        
        ClusterSimilarityMatrix.__Calculcate_Matrix__(matrix, \
            cluster_index_map, clusters, item_count, calculate)
                
        return ClusterSimilarityMatrix(matrix, cluster_index_map, item_count)
    
    @staticmethod
    def BuildFrom(available_clusters: List[Cluster], matrix: ClusterSimilarityMatrix) -> ClusterSimilarityMatrix:
        
        # formula_calculator: Callable[[Cluster, Cluster], float] = lambda c1, c2: \
        #     matrix[c1, c2] if matrix.hasCluster(c1,c2) \
        #     else formula.cluster_simularity(c1, c2, len(c1.intersection(c2)), matrix.__item_count__)
        item_count = matrix.__item_count__
        def cache_formula(c1: Cluster, c2: Cluster) -> float:
            if matrix.hasCluster(c1, c2):
                return matrix[c1, c2]
            return ClusterSimilarityMatrix.cluster_simularity(c1, c2, len(c1.intersection(c2)), item_count)
    
        size = len(available_clusters)
        cluster_index_map = {available_clusters[index]: index for index in range(len(available_clusters)) }
        new_matrix = np.full(shape = (size, size), fill_value=np.NINF, dtype=np.float32)
        
        new_matrix = ClusterSimilarityMatrix.__Calculcate_Matrix__(\
            new_matrix, cluster_index_map, available_clusters, item_count, cache_formula)
        return ClusterSimilarityMatrix(new_matrix, cluster_index_map, item_count)

    @staticmethod
    def __Calculcate_Matrix__(matrix: np.matrix, index_map: Dict[Cluster, int], \
        clusters: List[Cluster], item_count: int, calculation: Callable[[Cluster, Cluster], float]) -> np.matrix:
        for i in tqdm(range(len(clusters))):
            cluster1 = clusters[i]
            matrix_index_1 = index_map[cluster1]
            for j in range(i, len(clusters)):
                cluster2 = clusters[j]
                matrix_index_2 = index_map[cluster2]
                if i == j or cluster1 is cluster2 or cluster1.SamePartitionAs(cluster2):
                    # matrix[matrix_index_1, matrix_index_2] = np.NINF
                    continue
                
                value = calculation(cluster1, cluster2)
                
                matrix[matrix_index_1, matrix_index_2] = value
                matrix[matrix_index_2, matrix_index_1] = value
        return matrix
    
    
    @staticmethod
    def cluster_simularity(cluster1: Cluster, cluster2: Cluster, total_elements: int) -> float:    
        if cluster1.SamePartitionAs(cluster2): return np.NINF
        intersection_len = len(cluster1.intersection(cluster2))
        counter = intersection_len - ((len(cluster1) * len(cluster2)) / total_elements)
        divisor = sqrt(len(cluster1) * len(cluster2) * (1 - (len(cluster1) / total_elements)) * (1 - (len(cluster2) / total_elements)))
        return counter / divisor if divisor != 0 else 0
    
    def __len__(self) -> int:
        return len(self.index_map)

        
    def __getitem__(self, tuple: Tuple[Cluster, Cluster]) -> int or np.NINF:
        c1, c2 = tuple
        return self.getEntry(c1, c2)
    
    def __setitem__(self, object: Cluster) -> None:
        raise Exception("SimilarityMatrix is not allowed to be updated")
    
    
    def __expand__(self, cluster: Cluster, values: Dict[int, float or np.NINF]) -> None:
        Assert.assert_key_exists(cluster, self.index_map)
        new_index = self.index_map[cluster]
        size = len(values)
        
        new_row = np.empty_like(np.NINF, shape=(1, size-1), dtype=np.float32)
        new_column = np.empty_like(np.NINF, shape=(size, 1), dtype=np.float32)
        
        self.matrix = np.append(self.matrix, new_row,    axis=0 )       
        self.matrix = np.append(self.matrix, new_column, axis=1 )
        
        for index, value in values.items():
            self.matrix[new_index, index] = value
            self.matrix[index, new_index] = value
        return
    
    def ToDict(self, for_cluster: Cluster) -> Dict[Cluster, int or np.NINF]:
        Assert.assert_key_exists(for_cluster, self.index_map)
        index = self.index_map[for_cluster]
                
        return { cluster: self.getByIndex(index, other_index) for cluster, other_index in self.index_map.items()  }     
    
    
    def getEntry(self, cluster1 : Cluster, cluster2: Cluster) -> int or np.NINF:
        id1, id2 = self.index_map[cluster1], self.index_map[cluster2]
        return self.getByIndex(id1, id2)
    
    def updateCluster(self, cluster: Cluster) -> None:
        Assert.assert_not_none(cluster)
        Assert.assert_key_not_exists(cluster, self.index_map)
        index1 = self.index_map[cluster]
        for cluster2, index2 in self.index_map.items():
            if np.isneginf(self.matrix[index1, index2]):
                continue
            value = ClusterSimilarityMatrix.cluster_simularity(cluster, cluster2, cluster.intersection(cluster2), len(self.__item_count__))
            self.matrix[index1, index2] = value
            self.matrix[index2, index1] = value
    
    def add_merged_cluster(self, cluster: Cluster) -> None:
        Assert.assert_key_not_exists(cluster, self.index_map)
        
    
    def getByIndex(self, index1: int, index2: int) -> int or np.NINF:
        if index1 == index2:
            return np.NINF
        return self.matrix[index1, index2]
    
    def hasCluster(self, cluster1: Cluster, cluster2: Cluster) -> bool:
        return cluster1 in self.index_map and cluster2 in self.index_map
    
    def max_simularity_np(self) -> float:
        return self.matrix.max()
    
    def max_similarity(self, clusters: List[Cluster] = None) -> float:
        if clusters is None:
            return self.matrix.max()
        
        max_value = np.NINF
        for cl_index1 in range(len(clusters)):
            index1 = self.index_map[clusters[cl_index1]]
            for cl_index2 in range(cl_index1, len(clusters)):
                index2 = self.index_map[clusters[cl_index2]]
                value = self.getByIndex(index1, index2)
                if np.isneginf(value):
                    continue
                elif np.isneginf(max_value):
                    max_value = value 
                else:
                    max_value = max(value, max_value)
        return max_value
            
        