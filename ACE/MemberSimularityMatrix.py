from __future__ import annotations
import numpy as np
from tqdm import tqdm
from itertools import product
from Cluster import Cluster, PartitionSet
from typing import Iterable, Tuple, Dict, List
import Assertions as Assert 
import scipy.sparse as sp  
from ClusterSimilarityMatrix import DoubleSparseDictHashMatrix, SortKeysByHash, SparseDictHashMatrix, SparseTupleHashMatrix

def average(ite :Iterable) -> float:
    return sum(ite) / len(ite)


def Build_simularity_matrix(clusters: List[Cluster], gamma: PartitionSet) -> sp.dok_matrix:
    
    items = list(gamma.get_all_elements().keys())
    shape = (items, clusters)
    matrix = sp.dok_matrix(shape, dtype=np.float16)
    for cluster in clusters:
        membership_map = cluster.calc_all_membersimularity(max_member_value=len(gamma))
        for item, simularity in membership_map.items():
            matrix[item, cluster] = simularity
    return matrix


class MemberSimularityMatrix(DoubleSparseDictHashMatrix[object, Cluster, float]):
    def __init__(self, item_lst: List[object], cluster_lst: List[Cluster]) -> None:
        self.all_clusters = cluster_lst
        self.all_items = item_lst
        super().__init__(default_value=0, sparse_value_predicate=lambda x: x == 0)         
    
    @staticmethod
    def IndependentBuild(clusters: List[Cluster], gamma: PartitionSet):
        cluster_lst = list(clusters)
        item_lst = list(gamma.get_all_elements().keys()) 
        
        matrix = MemberSimularityMatrix(item_lst, cluster_lst)
        # matrix = sp.csc_matrix(shape, dtype=np.float32)
        
        for cluster in cluster_lst:
            membership_map = cluster.calc_all_membersimularity(max_member_value=len(gamma))
            for item, simularity in membership_map.items():
                matrix[item, cluster] = simularity
        return matrix
    
    def set(self, __k: Tuple[object, Cluster], __v: float):
        if __v == 0.0 and self.has_tuple(__k): 
            self.pop(__k)
            return
        if(__v > 1 or __v < 0):
            raise ValueError(f"similarity value '{__v}' is outside range 0-1")
        super().set(__k, __v)
            
    def Item_sum(self, item) -> float:
        return sum(self.get_row(item).values())
    
    def item_max(self, item) -> float:
        return max(self.get_row(item).values())
    
    def item_argMax(self, item) -> Tuple[Cluster, float]:
        arg_max, max_value = None, -1 
        
        row = self.get_row(item)
        if row is None or len(row) == 0:
            return (None, 0)
        
        for cluster, sim in row.items():
            if sim >= max_value:
                arg_max, max_value = cluster, sim
        if arg_max is None:
            raise KeyError("Could not find index in argmax")
        return (arg_max, max_value)
    
    def cluster_sum(self, cluster: Cluster) -> float:
        return sum(self.get_column(cluster).values())
    
    def cluster_max(self, cluster: Cluster) -> float:
        return sum(self.get_column(cluster).values())
        
    def cluster_mean(self, cluster: Cluster) -> float:
        return self.cluster_sum(cluster) / len(self.get_column(cluster))
        
    def assign_item_to(self, cluster: Cluster, item: object) -> None:
        for cluster2 in self.get_row(item).keys():
            if cluster is cluster2:
                continue
            self.pop_entry(item, cluster)
    
    def get_row(self, item: object) -> Dict[Cluster, float]:
        if item not in self.__internal_row__: return {}
        return dict(self.__internal_row__[item])
    
    def get(self, __k: Tuple[object, Cluster]) -> float:
        if self.has_tuple(__k):
            return super().get(__k)
        return 0.0


class MemberMatrix:
    def __init__(self, cluster_index_map: set[Cluster], item_index_map: set[object]) -> None:
        self.items = item_index_map
        self.clusters = cluster_index_map
        self.common_neighbor_cache: Dict[object, List[object]] = None
        
        
    @staticmethod
    def build(cluster_lst: List[Cluster], item_lst: List[object]) -> MemberMatrix:
        item_index_map = set(item_lst)
        cluster_index_map = set(cluster_lst)
        
        return MemberMatrix(cluster_index_map, item_index_map)
    
    
    # def get_all_common_items(self, totally_uncertain_items: set[object]) -> SparseDictHashMatrix[object, set[object]]:
    #     all_common_neighbors = SparseDictHashMatrix[object, set](default_value=[])

    #     for cluster in tqdm(self.clusters):
    #         intersection = totally_uncertain_items.intersection(cluster)
    #         for item1, item2 in [ (x, y) for x in intersection for y in cluster if x is not y]:
    #             if all_common_neighbors.has_entry(item1, item2) is False: all_common_neighbors[item1, item2] = set() 
    #             for item in [x for x in cluster if x is not item1 and x is not item2]:
    #                 all_common_neighbors[item1, item2].add(item)
    #                 all_common_neighbors[item2, item1].add(item)

    #     return all_common_neighbors
    
    def calculate_coassosiation_matrix(self, totally_uncertain_items: set[object], gamma: PartitionSet) -> SparseDictHashMatrix[object, Tuple[float, int]]:
        all_common_neighbors = SparseDictHashMatrix[object, Tuple[float, int]](SortKeysByHash, default_value=(0,0))
        dct = {item: gamma.calc_all_coassosiation(item) for item in totally_uncertain_items}

        def get(item1, item2) -> float:
            if item1 in dct and item2 in dct[item1]:
                return dct[item1][item2]
            if item2 in dct and item1 in dct[item2]:
                return dct[item2][item1]
            return 0.0

        for cluster in tqdm(self.clusters):
            intersection = totally_uncertain_items.intersection(cluster)
            for item1, item2 in [ (x, y) for x in intersection for y in cluster if x is not y and y not in intersection]:
                if all_common_neighbors.has_entry(item1, item2) is False: all_common_neighbors[item1, item2] = (0, 0)
                for common_item in [x for x in cluster if x is not item1 and x is not item2]:
                    value, count = all_common_neighbors[item1, item2]
                    all_common_neighbors[item1, item2] = (value + get(item1, common_item) + get(item2, common_item), count+1)
                    
        
        return all_common_neighbors
    
    
    def average_common_neighbors(self, item: object, cluster: Cluster, common_coassosiation_matrix: SparseDictHashMatrix[object, Tuple[float, int]] ) -> float:
        if len(cluster) <= 0:
            return 0
        if len(cluster) == 1 and item in cluster:
            return 0
        
        sum_value = sum([sum / (2*count) for sum, count in [common_coassosiation_matrix[item, item2] \
            for item2 in cluster if item != item2 ] if count != 0 ])
        
        return sum_value / len(cluster)
    
    # def initialize_cocache(self, totally_uncertain_items: set[object], coassociation_matrix: CoAssosiationMatrix) -> None:
    #     all_common_neighbors: Dict[object, float] = {}
        
    #     def add(item, cluster):
    #         if item not in all_common_neighbors: all_common_neighbors[item] = []
    #         all_common_neighbors[item].append(cluster)
        
    #     for cluster in tqdm(self.clusters):
    #         for item in totally_uncertain_items.intersection(cluster):
    #             add(item, cluster)

    #     self.common_neighbor_cache = all_common_neighbors
    
    # def initialize_cache(self, totally_uncertain_items: set[object]) -> None:
    #     all_common_neighbors: Dict[object, List[Cluster]] = {}
        
    #     def add(item, cluster):
    #         if item not in all_common_neighbors: all_common_neighbors[item] = []
    #         all_common_neighbors[item].append(cluster)
        
    #     for cluster in tqdm(self.clusters):
    #         for item in totally_uncertain_items.intersection(cluster):
    #             add(item, cluster)

    #     self.common_neighbor_cache = all_common_neighbors
    
    # def get_common_items(self, item1: object, item2: object) -> set[object]:
    #     Assert.assert_index_exists(item1, self.items)
    #     Assert.assert_index_exists(item2, self.items)
    #     # Assert.assert_not_none(self.common_neighbor_cache)
    #     # Assert.assert_key_exists(item1, self.common_neighbor_cache)
    #     Assert.assert_not_equal(item1, item2)
        
    #     common_neighbors = set()
    #     for cluster in self.clusters:
    #         if item1 in cluster and item2 in cluster: 
    #             common_neighbors.update(cluster.__iter__())
    #     return common_neighbors
    
    # def common_neighbors(self, coassociation_matrix: CoAssosiationMatrix, item1: object, item2: object,\
    #     common_items_matrix: SparseDictHashMatrix[object, set[object]] = None) -> float:
    #     #

    #     common_items = self.get_common_items(item1, item2) 
    #     if common_items is None or len(common_items) == 0:
    #         return 0
    #     sum_value = sum([coassociation_matrix[item1, c_item] + coassociation_matrix[c_item, item2] for c_item in common_items])
    #     div = 2* len(common_items)
    #     return sum_value / div if div != 0 else 0

    # def average_common_neighbors(self, coassociation_matrix: CoAssosiationMatrix, item: object, cluster: Cluster,\
    #     common_items_matrix: SparseDictHashMatrix[object, set[object]] = None) -> float:
    #     if len(cluster) <= 0:
    #         return 0
    #     if len(cluster) == 1 and item in cluster:
    #         return 0
    #     sumvalue = sum([self.common_neighbors(coassociation_matrix, item, item2, common_items_matrix) for item2 in cluster if item != item2 ])
    #     return sumvalue / len(cluster)
        

class CoAssosiationMatrix(SparseDictHashMatrix):
    def __init__(self) -> None:
        super().__init__(SortKeysByHash)
        # self.matrix = matrix
        # self.index_map = index_map
        
    
    @staticmethod
    def build_old(gamma: PartitionSet) -> CoAssosiationMatrix:
        item_lst = list((gamma.get_all_elements()).keys())
        index_map = {item_lst[i]: i for i in range(len(item_lst))}
        partition_count = len(gamma)
        shape = (len(item_lst), len(item_lst))
        # matrix :np.matrix = np.full(shape=shape, fill_value=0, dtype=np.float32)
        matrix = sp.lil_matrix(shape, dtype=np.float16)
        
        for item1 in tqdm(item_lst):
            matrix_index1 = index_map[item1]
            item_coasssiation =  gamma.calc_all_coassosiation(item1)
            
            for item2, value in item_coasssiation.items():
                matrix_index2 = index_map[item2]
                    
                matrix[matrix_index1, matrix_index2] = value
                matrix[matrix_index2, matrix_index1] = value
        
        return CoAssosiationMatrix(matrix, index_map, partition_count)
    
    @staticmethod
    def build(gamma: PartitionSet) -> CoAssosiationMatrix:
        item_lst = list((gamma.get_all_elements()).keys())
        matrix = CoAssosiationMatrix()
        
        for item1 in tqdm(item_lst):
            item_coasssiation =  gamma.calc_all_coassosiation(item1)
            for item2, value in item_coasssiation.items():
                matrix[item1, item2] = value
        
        return matrix
    
    def FindAssociatedItem(self, item: object) -> object or None:
        row = self.get_row(item)
        max_value, max_item = 0, None
        for other_item, value in row.items():
            if item is other_item: continue
            if value > max_value:
                max_value = value
                max_item = other_item
        return max_item
    
    def get(self, __k: Tuple[object, object]) -> float:
        if self.has_tuple(__k):
            return super().get(__k)
        return 0.0
        
# def partial_build_coassosiation_row(gamma: PartitionSet) -> Tuple[int, Dict[object, float]]:
#     pass