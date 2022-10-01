from __future__ import annotations
import itertools
from typing import Dict, List, Generic, TypeVar, Tuple, Iterable, Iterator
import Assertions as Assert
import sys
import numpy as np

sys.setrecursionlimit(10**6)

T = TypeVar("T")

class Cluster(Generic[T]):
    def __init__(self, partition_id = None):
        self.__partition_id__ = partition_id
        #maps each item in cluster to membership value
        self.__membership__: Dict[T, int] = dict()
        self.__total_size__: int = 0

    @staticmethod
    def merge(cluster1: Cluster[T], cluster2: Cluster[T]) -> Cluster[T]:
        Assert.assert_not_none(cluster1)
        Assert.assert_not_none(cluster2)
        cluster = Cluster(partition_id=None)
        
        for item, membership in itertools.chain(cluster1.__membership__.items(),\
                cluster2.__membership__.items()):
            cluster.__add_member__(item, membership)
        return cluster
    
    def __contains__(self, __o: T) -> bool:
        return __o in self.__membership__
    
    def __len__(self) -> int:
        return len(self.__membership__)
    
    def __iter__(self) -> Iterator[T]:
        return self.__membership__.__iter__()
    
    def intersection(self, other: Cluster) -> List[T]:
        #return list(set(self.__membership__.keys()) & set(other.__membership__.keys()))
        a, b = ((self, other) if len(self) < len(other) else (other, self))
        return set(a).intersection(b)
        # return [item for item in b if item in a]     
        
    def intersection_len(self, other: Cluster) -> int:
        return len(self.__membership__.keys() & other.__membership__.keys() )
        
    def union(self, other: Cluster) -> set[T]:
        return set(self).union(other)

    def append(self, __object: T) -> None:
        return self.add(__object)
    
    def add(self, __object: T) -> None:
        if __object in self.__membership__:
            return
            # raise Exception('Cluster already has item ' + str(__object))
        self.__membership__[__object] = 1

    def __add_member__(self, __object: T, membership: int) -> None:
        self.__membership__[__object] = self.__membership__[__object] + membership\
            if __object in self.__membership__ else membership

    def remove(self, __value: T) -> bool:
        if __value in self:
            self.__membership__.pop(__value)
            return True
        return False

    def descriminatory_union(self, other: Cluster[T]) -> List[T]:
        result = []
        for item in other:
            if item not in self:
                result.append(item)
        for item in self:
            if item not in other:
                result.append(item)
        return result

    def SamePartitionAs(self, other: Cluster) -> bool:
        if self.__partition_id__ is None or other.__partition_id__ is None:
            return False
        return self.__partition_id__ == other.__partition_id__

    def to_hash_vector(self) -> List[int]:
        return np.array([id(x) for x in self])

    def calc_all_membersimularity(self, max_member_value: int) -> Dict[T, float]:
        return { item: membersip / max_member_value for item, membersip in self.__membership__.items()}
    
    def calc_all_membership(self) -> Dict[T, int]:
        return self.__membership__

    def member_simularity(self, item: T, max_member_value: int) -> float:
        return self.membership(item) / max_member_value

    def mean_member_simularity(self, max_member_value: int) -> float:
        return self.sum_member_simularity(max_member_value) / len(self)
    
    def sum_member_simularity(self, max_member_value: int) -> float:
        return sum([self.member_simularity(item, max_member_value) for item in self])
    
    def max_member_simularity(self, max_member_value: int) -> float:
        return self.max_membership() / max_member_value
    
    def argmax_member_simularity(self, max_member_value: int) -> Tuple[T, float]:
        item, membership = self.argmax_membership()
        return (item, membership / max_member_value)
    
    def membership(self, item: T) -> int:
        return self.__membership__[item] if item in self.__membership__ else 0

    def mean_membership(self) -> float:
        return sum(self.__membership__.values()) / len(self)
    
    def sum_membership(self) -> float:
        return sum(self.__membership__.values())
    
    def max_membership(self) -> int:
        return max(self.__membership__.values())

    def argmax_membership(self) -> Tuple[T, int]:
        return max(self.__membership__.items(), key=lambda x: x[1])

    def __hash__(self) -> int:
        return id(self)
    
# Do not use normal set item of this dict, use add instead
class Partition(Dict[str, Cluster[T]], Generic[T]):

    def __init__(self, data: List[T] = None) -> None:
        self.__data__ = set(data) if data is not None else set()
        self.__bypass_assert__ = False if data is not None else True

    def add(self, cluster_name: str, item: T) -> None:
        if not self.__bypass_assert__: Assert.assert_partition_content(self.__data__, item)
        # self.__assert_item_not_in_other_cluster__(item)

        if cluster_name not in self:
            super().__setitem__(cluster_name, Cluster(id(self)))

        self[cluster_name].add(item)

    # def __setitem__(self, __k, __v) -> None:
        # raise Exception("DO NOT USE THIS, use '.add()' instead")
    def remove(self, cluster: Cluster[T]) -> None:
        del_name = None
        for name, cls2 in self.items(): 
            if cls2 is cluster:
                del_name = name
                break
        if del_name is not None:
            self.pop(del_name)
        
    def remove_lst(self, cluster_lst: List[Cluster[T]]) -> None:
        del_lst = []
        cluster_set = set(cluster_lst)
        for name, cls2 in self.items(): 
            if cls2 is cluster_set:
                del_lst.append(name)
        for del_name in del_lst:
            self.pop(del_name)
            

    def __assert_item_not_in_other_cluster__(self, item: T) -> None:
        for key, value in self.items():
            if item in value:
                raise Exception(f"Item {str(item)} is already in another cluster in this partition")

    def element_count(self) -> int:
        return len(self.__data__)
    
    def IsInSameCluster(self, item1: T, item2: T) -> bool:
        Assert.assert_item_in_list(self.__data__, item1)
        Assert.assert_item_in_list(self.__data__, item2)
        for cluster in self.values():
            if item1 in cluster and item2 in cluster:
                return True
        return False
    
    def get_neighboring_items(self, item: T) -> List[T]:
        Assert.assert_item_in_list(self.__data__, item)
        for cluster in self.values():
            if item in cluster:
                return list(cluster)
        return list([])
        

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

    def ClusterToPartitionMap(self) -> Dict[Cluster, int]:
        dct_info = {}
        for partition_idx in range(len(self)):
            for key, value in self[partition_idx].items():
                dct_info[value] = partition_idx
        return dct_info

    def calc_all_coassosiation(self, item: object) -> Dict[object, float]:
        Assert.assert_item_in_list(self.__dataset__, item)
        item_values = {}
        for partition in self:
            common = partition.get_neighboring_items(item)
            for d_item in common:
                item_values[d_item] = item_values[d_item] + 1 if d_item in item_values else 1
                    
        
        return { item: value / len(self) for item, value in item_values.items()  }

    def calc_coassosiation(self, item1: object, item2: object):
        Assert.assert_item_in_list(self.__dataset__, item1)
        Assert.assert_item_in_list(self.__dataset__, item2)
        
        return sum([1 if p.IsInSameCluster(item1, item2) else 0 for p in self]) / len(self)

    def get_all_clusters(self) -> List[Cluster[T]]:
        result = []
        for partition in self:
            for key, value in partition.items():
                result.append(value)
        
        return result

    def count_elements(self) -> int:
        return len(self.__dataset__)

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
    
    def get_all_items(self) -> List[T]:
        return list(self.__dataset__)

    def mean_cluster_in_partition(self) -> float:
        return sum( (len(partition) for partition in self) ) / len(self)

    def maximal_partition_clusters(self) -> int:
        return max([len(partition) for partition in self])
        
    
    
    
