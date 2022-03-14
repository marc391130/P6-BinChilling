from __future__ import annotations
from typing import Dict, List, Generic, TypeVar, Tuple
from math import sqrt
from typing_extensions import Self
import Assertions as Assert
import sys

sys.setrecursionlimit(10**6)

T = TypeVar("T")

class Contig:
    def __init__(self, name: str, vector: List[float]):
        self.name = name
        self.vector = vector


class Cluster(List[T]):
    def append(self, __object: T) -> None:
        if __object in self:
            raise Exception(f"Item {str(__object)} already in cluster")
        return super().append(__object)

    def intersection(self, other: Cluster[T]) -> Cluster[T]:
        result = []
        for item in other:
            if item in self:
                result.append(item)
        return result

    def contains(self, item: T) -> bool:
        return item in self

    def merge(self, other: Cluster[T]) -> Cluster[T]:
        result = Cluster()
        for x in self:
            result.append(x)
        for x in other:
            if x not in result:
                result.append(x)
        return result

    def __hash__(self) -> int:
        return id(self)

    @staticmethod
    def membership_similarity(cluster_lst: List[Cluster[T]], object: T) -> int:
        result = 0
        for cluster in cluster_lst:
            if cluster.contains(object):
                result += 1
        return result
    
# Do not use normal set item of this dict, use add instead
class Partition(Dict[str, Cluster[T]], Generic[T]):

    def __init__(self, data: List[T]) -> None:
        self.__data__ = data

    def add(self, cluster_name: str, item: T):
        Assert.assert_partition_content(self.__data__, item)
        self.__assert_item_not_in_other_cluster__(item)

        if cluster_name not in self:
            super().__setitem__(cluster_name, Cluster())

        self[cluster_name].append(item)

    def __setitem__(self, __k, __v) -> None:
        raise Exception("DO NOT USE THIS, use '.add()' instead")

    def __assert_item_not_in_other_cluster__(self, item: T) -> None:
        for key, value in self.items():
            if item in value:
                raise Exception(f"Item {str(item)} is already in another cluster in this partition")

    def element_count(self) -> int:
        return len(self.__data__)

class PartitionSet(List[Partition[T]]):
    def __init__(self):
        self.total_elements = 0

    def append(self, __object: Partition[T]) -> None:
        Assert.assert_partion_set_content(self)
        self.total_elements = max(self.total_elements, __object.element_count())
        return super().append(__object)

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


    def __similarity_measure_cluster__(self, cluster1: Cluster, cluster2: Cluster) -> float:
        cluster_intersection = cluster1.intersection(cluster2)

        counter = len(cluster_intersection) - ((len(cluster1) * len(cluster2)) / self.total_elements)
        divisor = sqrt(len(cluster1) * len(cluster2) * (1 - (len(cluster1) / self.total_elements)) * (1 - (len(cluster2) / self.total_elements)))
        return counter / divisor

    def get_all_clusters(self) -> List[Cluster[T]]:
        result = []
        for partition in self:
            for key, value in partition.items():
                result.append(value)
        return result

    def get_all_elements(self) -> List[T]:
        result = {}
        i = 0
        for item in self[0].__data__:
            result[item] = i if item not in result else result[item]
            i += 1
        return result

    def mean_cluster(self) -> float:
        return len(self.get_all_clusters()) / len(self)
    
    
