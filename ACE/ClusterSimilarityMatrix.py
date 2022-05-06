from __future__ import annotations
import itertools
from collections import Counter
import numpy as np
from tqdm import tqdm
from BinEvaluator import BinEvaluator
from Cluster import Cluster, PartitionSet
from multiprocessing import Pool, cpu_count
from typing import Iterable, Iterator, MutableMapping, Tuple, Dict, List, Callable, TypeVar, Generic
import Assertions as Assert 
from math import sqrt
import os
from Domain import ContigData
from SparseMatrix_implementations import SparseTupleHashMatrix, SparseDictHashMatrix, SortKeysByHash

TK = TypeVar("TK")
TV = TypeVar("TV")

UNIQUE_GENES = 104 

def cluster_simularity(cluster1: Cluster, cluster2: Cluster, total_elements: int) -> float:
    use_old = True 
    return cluster_simularity_old(cluster1, cluster2, total_elements) \
        if use_old else\
            cluster_simularity_new(cluster1, cluster2, total_elements) 
        

def cluster_simularity_old(cluster1: Cluster, cluster2: Cluster, total_elements: int) -> float:    
    if cluster1.SamePartitionAs(cluster2): return np.NINF
    intersection_len, len1, len2 = len(cluster1.intersection(cluster2)), len(cluster1), len(cluster2)
    counter = intersection_len - ((len1 * len2) / total_elements)
    divisor = sqrt(len1 * len2 * (1 - (len1 / total_elements)) * (1 - (len2 / total_elements)))
    return min(counter / divisor, 1.0) if divisor != 0 else 0.0


def cluster_simularity_new(cluster1: Cluster[ContigData], cluster2: Cluster[ContigData], total_elements: int) -> float:
    simularity = cluster_simularity_old(cluster1, cluster2, total_elements)
    
    com1, con1 = calc_con_fast(cluster1)
    com2, con2 = calc_con_fast(cluster2)
    comm, conm = calc_con_2(cluster1, cluster2)
    bonus = comm - max(com1, com2)
    penalty = (conm - min(con1, con2))
    # if penalty > 0: print(penalty)
    # elif penalty < 0: print(conm, ' - ', max(con1, con2), ' = ', penalty)
    return simularity + 0.0*(bonus - penalty)

def calc_con_fast(cluster: List[ContigData]) -> Tuple[float, float]:
    seen, dup = set(), set()
    for item in cluster:
        for scg in item.SCG_genes:
            if scg in seen:
                dup.add(scg)
            seen.add(scg)
    return (len(dup) / len(seen), len(seen) / UNIQUE_GENES) if len(seen) != 0 else (0.0, 0.0)

def calc_con_2(cluster1: List[ContigData], cluster2: List[ContigData]) -> float:
    seen_contigs, seen_scg, dup_scg = set(), set(), set()
    for item in itertools.chain(cluster1, cluster2):
        if item in seen_contigs: continue
        seen_contigs.add(item)
        for scg in item.SCG_genes:
            if scg in seen_scg:
                dup_scg.add(scg)
            seen_scg.add(scg)
    return (len(dup_scg) / len(seen_scg), len(seen_scg) / UNIQUE_GENES) if len(seen_scg) != 0 else (0.0, 0.0)

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
    def build_multithread_old(cluster_lst: List[Cluster], total_item_count: int, a1_min: float,\
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
    def build_multithread(cluster_lst: List[Cluster], total_item_count: int, a1_min: float,\
        processes: int, chunksize: int = 75) -> SparseClustserSimularity:
        if processes <= 1: return SparseClustserSimularity.build(cluster_lst, total_item_count, a1_min)
        matrix_dct = SparseTupleHashMatrix(SortKeysByHash)
        global shared_cluster_lst_memory 
        shared_cluster_lst_memory = cluster_lst
                
        parameters = [(i, total_item_count, a1_min) for i in range(len(cluster_lst))]
        # result : List[List[Tuple[Cluster, Cluster, float]]] = None
        with Pool(min(cpu_count(), processes)) as p:
            for chunk_result in tqdm(p.imap_unordered(partial_build_similarity_row_shared_mem, parameters, chunksize=chunksize), total=len(cluster_lst)):
                for i1, i2, sim in chunk_result:
                    matrix_dct[ cluster_lst[i1], cluster_lst[i2] ] = sim
        
        del shared_cluster_lst_memory
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

shared_cluster_lst_memory : List[Cluster] or None = None 
def partial_build_similarity_row_shared_mem(tup: Tuple[int, int, int]) -> List[Tuple[int, int, float]]:
    i1, total_count, a1_min = tup
    global shared_cluster_lst_memory
    cluster_lst = shared_cluster_lst_memory
    lst, c1 = [], cluster_lst[i1]
    for i2 in range(i1+1, len(cluster_lst)):
        c2 = cluster_lst[i2]
        similarity = cluster_simularity(c1, c2, total_count)
        if similarity >= a1_min:
            lst.append( (i1, i2, similarity) )
    return lst