from __future__ import annotations
from typing import Dict, List, Generic, TypeVar, Tuple
from math import sqrt
import Constants as Constant

import Assertions as Assert
import sys

sys.setrecursionlimit(10**6)

T = TypeVar("T")

class Contig:
    def __init__(self, name: str):
        self.name = name

    def __hash__(self) -> int:
        return id(self)
    
    def __str__(self):
        return self.name

class Cluster(set, Generic[T]):
    def __init__(self, partition_id = None):
        self.__children_lst__ = []
        self.__partition_id__ = partition_id

    @staticmethod
    def merge(cluster1: Cluster[T], cluster2: Cluster[T]) -> Cluster[T]:
        Assert.assert_not_none(cluster1)
        Assert.assert_not_none(cluster2)
        result = Cluster()
        result.__children_lst__ = [cluster1, cluster2]
        for x in cluster1:
            result.append(x)
        for x in cluster2:
            if x not in result:
                result.append(x)
        return result
    
    def append(self, __object: T) -> None:
        if __object in self:
            raise Exception(f"Item {str(__object)} already in cluster")
        return super().add(__object)

    def intersection(self, other: Cluster[T]) -> List[T]:
        result = []
        for item in other:
            if item in self:
                result.append(item)
        return result

    def remove(self, __value: T) -> None:
        if __value not in self:
            return
        super().remove(__value)
        for child in self.__children_lst__:
            child.remove(__value)

    def contains(self, item: T) -> bool:
        return item in self
    
    def SamePartitionAs(self, other: Cluster) -> bool:
        if self.__partition_id__ is None or other.__partition_id__ is None:
            return False
        return self.__partition_id__ == other.__partition_id__

    # def merge(self, other: Cluster[T]) -> Cluster[T]:
    #     result = Cluster(self, other)
    #     for x in self:
    #         result.append(x)
    #     for x in other:
    #         if x not in result:
    #             result.append(x)
    #     return result
    


    def calc_membership(self, item: T) -> int:
        item_membership_values = self.calc_all_membership()
        return item_membership_values[item] if item in item_membership_values else 0
    
    def calc_all_membership(self) -> Dict[T, int]:
        if len(self.__children_lst__) == 0:
            return {item: 1 for item in self}
        
        result = {}
        for child in self.__children_lst__:
            child_dic = child.calc_all_membership()
            for item, value in child_dic.items():
                result[item] = value if item not in result else result[item] + value
                
        for item in self:
            if item not in result:
                result[item] = 1
        Assert.assert_equal(len(result), len(self))
        return result


    def __get_leaf_clusters__(self) -> List[Cluster]:
        result = []
        
        if len(self.__children_lst__) == 0:
            return [self]

        for child in self.__children_lst__:
            result += child.__get_leaf_clusters__()
        
        return result

    def __hash__(self) -> int:
        return id(self)

    # @staticmethod
    # def membership_similarity(cluster_lst: List[Cluster[T]], object: T) -> int:
    #     result = 0
    #     for cluster in cluster_lst:
    #         if cluster.contains(object):
    #             result += 1
    #     return result
    
# Do not use normal set item of this dict, use add instead
class Partition(Dict[str, Cluster[T]], Generic[T]):

    def __init__(self, data: List[T]) -> None:
        self.__data__ = data

    def add(self, cluster_name: str, item: T) -> None:
        Assert.assert_partition_content(self.__data__, item)
        self.__assert_item_not_in_other_cluster__(item)

        if cluster_name not in self:
            super().__setitem__(cluster_name, Cluster(id(self)))

        self[cluster_name].append(item)

    # def __setitem__(self, __k, __v) -> None:
        # raise Exception("DO NOT USE THIS, use '.add()' instead")

    def __assert_item_not_in_other_cluster__(self, item: T) -> None:
        for key, value in self.items():
            if item in value:
                raise Exception(f"Item {str(item)} is already in another cluster in this partition")

    def element_count(self) -> int:
        return len(self.__data__)

class PartitionSet(List[Partition[T]]):
    def __init__(self, data: List[T]):
        Assert.assert_list_nonempty(data)
        self.__dataset__ = data

    def append(self, __object: Partition[T]) -> None:
        Assert.assert_partion_set_content(self)
        self.total_elements = max(self.total_elements, __object.element_count())
        return super().append(__object)
    
    def create_partition(self) -> Partition[T]:
        partition = Partition(self.__dataset__)
        super().append(partition)
        return partition

    def similarity_measure(self, gamma_idx1: int, gamma_idx2: int, cluster1_name: str, cluster2_name: str) -> float: 
        return self.similarity_measure((gamma_idx1, cluster1_name), (gamma_idx2, cluster2_name))

    def similarity_measure(self, cluster1_idx: Tuple[int, str], cluster2_idx: Tuple[int, str]) -> float: 
        Assert.assert_index_exists(cluster1_idx[0], self)
        Assert.assert_index_exists(cluster2_idx[0], self)

        partition1 = self[cluster1_idx[0]]
        partition2 = self[cluster2_idx[0]]

        Assert.assert_key_exists(cluster1_idx[1], partition1)
        Assert.assert_key_exists(cluster2_idx[1], partition2)

        cluster1 = partition1[cluster1_idx[1]]
        cluster2 = partition2[cluster2_idx[1]]

        return self.__similarity_measure_cluster__(cluster1, cluster2)

    def ClusterToPartitionMap(self) -> Dict[Cluster, int]:
        dct_info = {}
        for partition_idx in range(len(self)):
            for key, value in self[partition_idx].items():
                dct_info[value] = partition_idx
        return dct_info


    def get_all_clusters(self) -> List[Cluster[T]]:
        result = []
        for partition in self:
            for key, value in partition.items():
                result.append(value)
        
        return result

    def __total_elements__(self) -> int:
        return len(self.__dataset__)

    # def get_all_elements(self) -> List[T]:
    #     return list(self.__dataset__).copy()
    
    def get_all_elements(self) -> Dict[T, int]:
        result = {}
        i = 0
        for item in self.__dataset__:
            result[item] = i if item not in result else result[item]
            i += 1
        return result

    def mean_cluster_in_partition(self) -> float:
        return len(self.get_all_clusters()) / len(self)

    def maximal_partition_clusters(self) -> int:
        return max([len(partition) for partition in self])
        
    
    
    
