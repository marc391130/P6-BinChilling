from __future__ import annotations
import itertools
import multiprocessing
import numpy as np
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from multiprocessing import Pool, get_context
from typing import Iterable, Iterator, MutableMapping, Tuple, Dict, List, Callable, TypeVar, Generic
import Assertions as Assert 
from math import sqrt

TK = TypeVar("TK")
TV = TypeVar("TV")

def cluster_simularity(cluster1: Cluster, cluster2: Cluster, total_elements: int) -> float:    
    if cluster1.SamePartitionAs(cluster2): return np.NINF
    intersection_len = len(cluster1.intersection(cluster2))
    counter = intersection_len - ((len(cluster1) * len(cluster2)) / total_elements)
    divisor = sqrt(len(cluster1) * len(cluster2) * (1 - (len(cluster1) / total_elements)) * (1 - (len(cluster2) / total_elements)))
    return counter / divisor if divisor != 0 else 0


def SortKeys(key1: TK, key2: TK) -> Tuple[TK, TK]:
    return (key1, key2) if key1.__hash__() <= key2.__hash__() else (key2, key1)

def SortTuple(__k: Tuple[TK, TK]) -> Tuple[TK, TK]:
    return SortKeys(__k[0], __k[1])

class HashIterator(Iterator[Tuple[TK, TV]]):
    def __init__(self, source: Dict[Tuple[TK, TK], TV], pivot: TK, control: Iterator[TK], second_iterator: Iterator[TK] = None) -> None:
        self.source = source
        self.pivot = pivot
        self.control = control if second_iterator is None else itertools.chain(control, second_iterator)
        
    def __iter__(self) -> Iterator[Tuple[TK, TV]]:
        return self
        
    def __next__(self) -> Tuple[TK, TV]:
        tupkey = SortKeys(self.pivot, self.control.__next__())
        return self.source[tupkey]

class SparseHashMatrix(MutableMapping[Tuple[TK, TK], TV]):
    def __init__(self) -> None:
        self.__internal__ : Dict[Tuple[TK, TK], TV] = dict()
            
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK) -> TV:
        return self.get( SortKeys(key1, key2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK]]:
        return self.__internal__.keys()
    
    def values(self) -> Iterable[TV]:
        return self.__internal__.values()
    
    def items(self) -> Iterable[Tuple[TK, TK], TV]:
        return self.__internal__.items()
    
    def get(self, __k: Tuple[TK, TK]) -> TV:
        return self.__internal__[SortTuple(__k)]
    
    def __getitem__(self, __k: Tuple[TK, TK]) -> TV:
        return self.get(__k)
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK], __v: TV) -> None:
        tup = SortTuple(__k)
        self.__internal__[tup] = __v
        
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK]) -> None:
        return self.pop(__v)
    
    def pop(self, __k: Tuple[TK, TK]) -> None:
        self.__internal__.pop(__k)
    
    def pop_entry(self, key1: TK, key2: TK) -> None:
        tup = SortKeys(key1, key2)
        Assert.assert_key_exists(tup, self.__internal__)
        self.__internal__.pop(tup)

    def pop_contain(self, key: TK) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] is key or tup[1] is key: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    def pop_set(self, keys: set[TK]) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] in keys or tup[1] in keys: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    #UTILITY FUNCTIONS
    def has_entry(self, key: Tuple[TK, TK]) -> bool:
        tup = SortTuple(key)
        return tup in self.__internal__
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        return False
    
    def __len__(self) -> int:
        return len(self.__internal__)
    
    def __iter__(self) -> Iterator[Tuple[TK, TK], TV]:
        return self.__internal__.__iter__()

    
class SparseClustserSimularity:
    def __init__(self, cluster_lst: List[Cluster], total_item_count: int, min_value: float, matrix: SparseHashMatrix[Cluster, float] = None) -> None:
        self.matrix = SparseHashMatrix[Cluster, float]() if matrix is None else matrix
        self.__all_clusters__ = set(cluster_lst)
        self.__non_similar_clusters__ = set(cluster_lst)
        self.total_item_count = total_item_count
        self.min_value = min_value

    @staticmethod
    def build(cluster_lst: List[Cluster], total_item_count: int, a1_min: float) -> SparseClustserSimularity:
        matrix = SparseClustserSimularity(cluster_lst, total_item_count, a1_min)
        cluster_len = len(cluster_lst)
        
        for c1 in tqdm(range(cluster_len)):
            cluster1 = cluster_lst[c1]
            for c2 in range(c1+1, cluster_len):
                cluster2 = cluster_lst[c2]
                simularity = cluster_simularity(cluster1, cluster2, total_item_count)
                if simularity > a1_min:
                    matrix.set_value(cluster1, cluster2, simularity)
        return matrix
    
    @staticmethod
    def build_multithread(cluster_lst: List[Cluster], total_item_count: int, a1_min: float, processes: int, chunksize: int = 75) -> SparseClustserSimularity:
        if processes <= 1: return SparseClustserSimularity.build(cluster_lst, total_item_count, a1_min)
        matrix_dct = SparseHashMatrix()
        parameters = [(i, cluster_lst, total_item_count, a1_min) for i in range(len(cluster_lst))]
        with Pool(processes) as p:
            result : List[List[Tuple[Cluster, Cluster, float]]] =\
                tqdm(p.imap(partial_build_similarity_row, parameters, chunksize=chunksize), total=len(parameters))
            
            for i1, i2, sim in itertools.chain.from_iterable(result):
                matrix_dct[ cluster_lst[i1], cluster_lst[i2] ] = sim
            
        return SparseClustserSimularity(cluster_lst, total_item_count, a1_min, matrix_dct)
    
    
    @staticmethod
    def hash_sort_cluster(a: Cluster, b: Cluster) -> Tuple[Cluster, Cluster]:
        return (a, b) if a.__hash__() <= b.__hash__() else (b, a)

    def get_available_clusters(self) -> List[Cluster]:
        return list(self.__all_clusters__)

    def remove_cluster(self, cluster: Cluster) -> None:
        Assert.assert_key_exists(cluster, self.__all_clusters__)
        self.__all_clusters__.remove(cluster)
        self.__non_similar_clusters__.discard(cluster)
        self.matrix.pop_contain(cluster)
    
    def remove_set_of_clusters(self, clusters: set[Cluster]) -> None:
        for cluster in clusters:
            Assert.assert_key_exists(cluster, self.__all_clusters__)
            self.__all_clusters__.remove(cluster)
            self.__non_similar_clusters__.discard(cluster)
        self.matrix.pop_set(clusters)
    
    def add_cluster(self, cluster: Cluster):
        Assert.assert_key_not_exists(cluster, self.__all_clusters__)
        self.__all_clusters__.add(cluster)
        self.__non_similar_clusters__.add(cluster)
        
    
    def set_value(self, cluster1: Cluster, cluster2: Cluster, value: float) -> bool:
        Assert.assert_key_exists(cluster1, self.__all_clusters__)
        Assert.assert_key_exists(cluster2, self.__all_clusters__)
        Assert.assert_not_equal(cluster1, cluster2)
        if value < self.min_value: return False
        self.__non_similar_clusters__.discard(cluster1)
        self.__non_similar_clusters__.discard(cluster2)
        self.matrix.set( (cluster1, cluster2), value )
        return True
    
    # def ToList(self) -> Dict[Cluster, Dict[Cluster, float]]:
    #     return [ (first_key, [(second_key, self.matrix.__internal__[first_key, second_key]) ] )\
    #         for first_key, second_key in self.matrix.__internal_firstkey__.items()]
    
    def get_entries(self) -> Tuple[ Tuple[Cluster, Cluster], float]:
        return self.matrix.__internal__.items()
    
    def shape(self) -> Tuple[int, int]:
        return (len(self.__all_clusters__), len(self.__all_clusters__))
        
    def __getitem__(self, tuple: Tuple[Cluster, Cluster]) -> float or None:
        c1, c2 = tuple
        return self.getEntry(c1, c2)
    
    def getEntry(self, cluster1: Cluster, cluster2: Cluster) -> float or None:
        return self.tryGetEntry(cluster1, cluster2, None)
    
    def tryGetEntry(self, cluster1: Cluster, cluster2: Cluster, default_value = None) -> float:
        if self.matrix.has_entry( (cluster1, cluster2) ):
            return self.matrix[cluster1, cluster2]
        return default_value
    
def partial_build_similarity_row(tup: Tuple[int, List[Cluster], int, int]) -> List[Tuple[int, int, float]]:
    i1, cluster_lst, total_count, a1_min = tup
    lst = []
    c1 = cluster_lst[i1]
    for i2 in range(i1+1, len(cluster_lst)):
        c2 = cluster_lst[i2]
        similarity = cluster_simularity(c1, c2, total_count)
        if similarity >= a1_min:
            lst.append( (i1, i2, similarity) )
    return lst