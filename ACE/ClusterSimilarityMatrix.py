from __future__ import annotations
import itertools
import multiprocessing
import numpy as np
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from multiprocessing import Pool, cpu_count
from typing import Iterable, Iterator, MutableMapping, Tuple, Dict, List, Callable, TypeVar, Generic
import Assertions as Assert 
from math import sqrt
import os

TK = TypeVar("TK")
TV = TypeVar("TV")

def cluster_simularity(cluster1: Cluster, cluster2: Cluster, total_elements: int) -> float:    
    if cluster1.SamePartitionAs(cluster2): return np.NINF
    intersection_len, len1, len2 = len(cluster1.intersection(cluster2)), len(cluster1), len(cluster2)
    counter = intersection_len - ((len1 * len2) / total_elements)
    divisor = sqrt(len1 * len2 * (1 - (len1 / total_elements)) * (1 - (len2 / total_elements)))
    return counter / divisor if divisor != 0 else 0


def SortKeysByHash(key1: TK, key2: TK) -> Tuple[TK, TK]:
    return (key1, key2) if key1.__hash__() <= key2.__hash__() else (key2, key1)


class HashIterator(Iterator[Tuple[TK, TV]]):
    def __init__(self, source: Dict[Tuple[TK, TK], TV], keysort : Callable[[TK, TK], Tuple[TK, TK]],\
        pivot: TK, control: Iterator[TK], second_iterator: Iterator[TK] = None) -> None:
        self.source = source
        self.pivot = pivot
        self.control = control if second_iterator is None else itertools.chain(control, second_iterator)
        self.keysort = keysort
        
    def __iter__(self) -> Iterator[Tuple[TK, TV]]:
        return self
        
    def __next__(self) -> Tuple[TK, TV]:
        tupkey = self.keysort(self.pivot, self.control.__next__())
        return self.source[tupkey]

class SparseDictHashMatrix(MutableMapping[Tuple[TK, TK], TV]):
    def __init__(self, keysort: Callable[[TK, TK], Tuple[TK, TK]] = None, default_value = None,\
        sparse_value_predicate: Callable[[TV], bool] = None) -> None:
        #sparses the value if the predicate returns true 
        
        self.__internal__ : Dict[TK, Dict[TK, TV]] = dict()
        self.keysort = keysort if keysort is not None else lambda x,y: (x,y)
        self.__default__ = default_value
        self.__sparse_filter__ = sparse_value_predicate if sparse_value_predicate is not None\
            else lambda x: x == self.__default__
        
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK) -> TV:
        k1, k2 = self.keysort(key1, key2)
        return self.get( (k1, k2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK]]:
        return self.__internal__.keys()
    
    def values(self) -> Iterable[TV]:
        return self.__internal__.values()
    
    def items(self) -> Iterable[Tuple[TK, Dict[TK, TV]]]:
        return self.__internal__.items()
    
    def get(self, __k: Tuple[TK, TK]) -> TV:
        k1, k2 = self.keysort(__k[0], __k[1])
        return self.__internal__.get(k1, {}).get(k2, self.__default__)
    
    def __getitem__(self, __k: Tuple[TK, TK]) -> TV:
        return self.get(__k)
    
    def get_row(self, __k: TK) -> Dict[TK, TV]:
        return self.__internal__.get(__k)
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK], __v: TV) -> None:
        if self.__sparse_filter__(__v): return
        
        k1, k2 = self.keysort(__k[0], __k[1])
        if k1 not in self.__internal__: self.__internal__[k1] = {}
        self.__internal__[k1][k2] = __v
    
    def set_dict(self, key: TK, dct: Dict[TK, TV]) -> None:
        for other_key, value in dct.items():
            self.set( self.keysort( key, other_key ), value)

    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK]) -> TV:
        k1, k2 = self.keysort(__v[0], __v[1])
        return self.__internal__[k1].pop(k2)
    
    def pop(self, __k: Tuple[TK, TK]) -> None:
        return self.__delitem__(__k)
    
    def pop_entry(self, key1: TK, key2: TK) -> None:
        return self.__delitem__( self.keysort(key1, key2) )

    def pop_contain(self, key: TK) -> None: #will always return none and never throw
        for key, inner_dct in self.__internal__:
            inner_dct.pop(key, None)
        self.__internal__.pop(key, None)
    
    #UTILITY FUNCTIONS
    def has_row_key(self, key: TK) -> bool:
        return key in self.__internal__
    
    def has_tuple(self, key: Tuple[TK, TK]) -> bool:
        k1, k2 = self.keysort(key[0], key[1])
        return self.has_entry(k1, k2)
    
    def has_entry(self, k1: TK, k2: TK) -> bool:
        k1, k2 = self.keysort(k1, k2)
        return k1 in self.__internal__ and k2 in self.__internal__[k1]
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        return False
    
    def __len__(self) -> int:
        return len(self.__internal__)
    
    def __iter__(self) -> Iterator[Tuple[TK, Dict[TK, TV]]]:
        return self.__internal__.__iter__()
    
TK2 = TypeVar('TK2')
    
class DoubleSparseDictHashMatrix(MutableMapping[Tuple[TK, TK2], TV]):
    def __init__(self, default_value = None, sparse_value_predicate: Callable[[TV], bool] = None) -> None:
        #sparses the value if the predicate returns true 
        
        self.__internal_row__ : Dict[TK, Dict[TK2, TV]] = dict()
        self.__internal_column__ : Dict[TK2, Dict[TK, TV]] = dict()
        self.__default__ = default_value
        self.__sparse_filter__ = sparse_value_predicate if sparse_value_predicate is not None\
            else lambda x: x == self.__default__
        
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK2) -> TV:
        return self.get( (key1, key2) )
    
    
    def items(self) -> Iterable[Tuple[TK, Dict[TK, TV]]]:
        return self.__internal_row__.items()
    
    def get(self, __k: Tuple[TK, TK2]) -> TV:
        k1, k2 = __k
        return self.__internal_row__.get(k1, {}).get(k2, self.__default__)
    
    def __getitem__(self, __k: Tuple[TK, TK2]) -> TV:
        return self.get(__k)
    
    def get_row(self, __k: TK) -> Dict[TK2, TV]:
        return self.__internal_row__.get(__k, {})
    
    def get_column(self, __k: TK2) -> Dict[TK, TV] or None:
        return self.__internal_column__.get(__k, {})
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK2], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK2], __v: TV) -> None:
        if self.__sparse_filter__(__v): return
        k1, k2 = __k
        if k1 not in self.__internal_row__: self.__internal_row__[k1] = {}
        self.__internal_row__[k1][k2] = __v
        if k2 not in self.__internal_column__: self.__internal_column__[k2] = {}
        self.__internal_column__[k2][k1] = __v
    
    def set_row(self, key: TK, dct: Dict[TK2, TV]) -> None:
        for other_key, value in dct.items():
            self.set( (key, other_key) , value)

    def set_column(self, key: TK2, dct: Dict[TK, TV]) -> None:
        for other_key, value in dct.items():
            self.set( ( other_key, key), value)
    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK2]) -> TV:
        k1, k2 = __v
        x1 = self.__internal_row__[k1].pop(k2, self.__default__)
        x2 = self.__internal_column__[k2].pop(k1, self.__default__)
        return x1
    
    def pop(self, __k: Tuple[TK, TK2]) -> None:
        return self.__delitem__(__k)
    
    def pop_entry(self, key1: TK, key2: TK2) -> None:
        return self.__delitem__( (key1, key2) )
    
    #UTILITY FUNCTIONS
    def __len__(self) -> Tuple[int, int]:
        return (len(self.__internal_row__), len(self.__internal_column__))
    
    def has_row_key(self, key: TK) -> bool:
        return key in self.__internal_row__
    
    def has_column_key(self, key: TK2) -> bool:
        return key in self.__internal_column__
    
    def has_tuple(self, key: Tuple[TK, TK2]) -> bool:
        k1, k2 = key
        return self.has_entry(k1, k2)
    
    def has_entry(self, k1: TK, k2: TK2) -> bool:
        return k1 in self.__internal_row__ and k2 in self.__internal_row__.get(k1, {})
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        return False
    
    def __iter__(self) -> Iterator[Tuple[TK, Dict[TK2, TV]]]:
        return self.__internal_row__.__iter__()
    
    
class SparseTupleHashMatrix(MutableMapping[Tuple[TK, TK], TV]):
    def __init__(self, keysort : Callable[[TK, TK], Tuple[TK, TK]] = None) -> None:
        self.__internal__ : Dict[Tuple[TK, TK], TV] = dict()
        self.keysort = keysort if keysort is not None else lambda x, y: (x, y) 
    
    def sortTuple(self, tup: Tuple[TK, TK]) -> Tuple[TK, TK]:
        return self.keysort(tup[0], tup[1])
    
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK) -> TV:
        return self.get( self.keysort(key1, key2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK]]:
        return self.__internal__.keys()
    
    def values(self) -> Iterable[TV]:
        return self.__internal__.values()
    
    def items(self) -> Iterable[Tuple[TK, TK], TV]:
        return self.__internal__.items()
    
    def get(self, __k: Tuple[TK, TK]) -> TV:
        return self.__internal__[self.sortTuple(__k)]
    
    def __getitem__(self, __k: Tuple[TK, TK]) -> TV:
        return self.get(__k)
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK], __v: TV) -> None:
        tup = self.sortTuple(__k)
        self.__internal__[tup] = __v
    
    def set_dict(self, key: TK, dct: Dict[TK, TV] ):
        for other_key, value in dct.items():
            self.set( (key, other_key), value )
    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK]) -> None:
        return self.pop(__v)
    
    def pop(self, __k: Tuple[TK, TK]) -> None:
        self.__internal__.pop(__k)
    
    def pop_entry(self, key1: TK, key2: TK) -> None:
        tup = self.keysort(key1, key2)
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
        tup = self.sortTuple(key)
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
    def __init__(self, cluster_lst: List[Cluster], total_item_count: int, min_value: float, matrix: SparseTupleHashMatrix[Cluster, float] = None) -> None:
        self.matrix = SparseTupleHashMatrix[Cluster, float](SortKeysByHash) if matrix is None else matrix
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
    def build_multithread(cluster_lst: List[Cluster], total_item_count: int, a1_min: float,\
        processes: int, chunksize: int = 75) -> SparseClustserSimularity:
        if processes <= 1: return SparseClustserSimularity.build(cluster_lst, total_item_count, a1_min)
        matrix_dct = SparseTupleHashMatrix(SortKeysByHash)
        
        parameters = [(i, cluster_lst, total_item_count, a1_min) for i in range(len(cluster_lst))]
        # result : List[List[Tuple[Cluster, Cluster, float]]] = None
        with Pool(min(cpu_count(), processes)) as p:
            for chunk_result in tqdm(p.imap(partial_build_similarity_row, parameters, chunksize=chunksize), total=len(cluster_lst)):
                for i1, i2, sim in chunk_result:
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
    
    def __len__(self):
        return len(self.__all_clusters__)
    
def partial_build_similarity_row(tup: Tuple[int, List[Cluster], int, int]) -> List[Tuple[int, int, float]]:
    i1, cluster_lst, total_count, a1_min = tup
    lst, c1 = [], cluster_lst[i1]
    for i2 in range(i1+1, len(cluster_lst)):
        c2 = cluster_lst[i2]
        similarity = cluster_simularity(c1, c2, total_count)
        if similarity >= a1_min:
            lst.append( (i1, i2, similarity) )
    return lst

def partial_build_similarity_row_double(tup: Tuple[int, List[Cluster], int, int]) -> List[Tuple[int, int, float]]:
    i1, cluster_lst, total_count, a1_min = tup
    lst, c1, total = [], cluster_lst[i1], len(cluster_lst)
    for i2 in range(i1+1, total-1, 2):
        c2, c3 = cluster_lst[i2], cluster_lst[i2+1]
        similarity1, similarity2 = cluster_simularity(c1, c2, total_count), cluster_simularity(c1, c3, total_count)
        if similarity1 >= a1_min: lst.append( (i1, i2, similarity1) )
        if similarity2 >= a1_min: lst.append( (i1, i2+1, similarity2) )
    if total % 2 == 1:
        c2 = cluster_lst[len(cluster_lst)-1]
        similarity = cluster_simularity(c1, c2, total_count)
        if similarity >= a1_min: lst.append( (i1, i2, similarity) )
    return lst