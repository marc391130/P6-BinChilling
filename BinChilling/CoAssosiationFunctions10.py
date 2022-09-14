import os
from sre_constants import OP_IGNORE
import tempfile
from typing import Callable, List, Dict, Tuple
from Cluster_matrices import MemberSimularityMatrix
from ClusterDomain import Cluster
import multiprocessing as mp
from SparseMatrix_implementations import SparseDictHashMatrix
from tqdm import tqdm
from Cluster_matrices import CoAssosiationMatrix
import numpy as np
from SharedDatastructures_implementations import SharedHashTable, UnevenNestedList
import sys
import gc
from math import floor, e
        
# (combined hash of two items) -> (co-assosiation value ( CO function ) for those two items)
shared_co_dct: SharedHashTable = None
shared_hash_lst: UnevenNestedList = None
memfile = tempfile.mktemp(suffix=".npz")

def read_shared_hash_item(index: int) -> List[int]:
    global shared_hash_lst
    return shared_hash_lst[index]

def recalculate_simularity_multiprocess(item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
    cluster_lst: List[Cluster], chunksize: int) -> MemberSimularityMatrix:
    
    global shared_co_dct
    if(shared_co_dct is None):
        raise Exception("shared memory co-matrix is not instantiated. call 'build_shared_memory_co_matrix' with reference to it, to instantiate it.")
    
    #Generate the reduced cluster lst map and cluster_hash_lst
    #reduced_cluster_lst_map is a map of object, to a list of indexes for clusters
        #this is used to not have to loop over all clusters, only the relevant ones.
    reduced_cluster_lst_map: Dict[object, List[int]] = {}
    cluster_hash_lst: List[List[int]] = np.empty(len(cluster_lst), dtype=object)
    for c_index, cluster in tqdm(enumerate(cluster_lst), total=len(cluster_lst), desc='Prep. and reducing param.  '):
        cluster_hash_lst[c_index] = np.array([hash(item) for item in cluster], dtype=np.int_)
        for item in cluster:
            reduced_cluster_lst_map[item] = \
                reduced_cluster_lst_map.get(item, []) + [c_index]
    
    parameters = ((i, hash(item), reduced_cluster_lst_map.get(item, [])) for i, item in tqdm(enumerate(item_lst), desc='Recalc. similarity matrix  ', total=len(item_lst)))
    
    try:
        with mp.Pool(mp.cpu_count(), initializer=init_mp, initargs=(shared_co_dct, cluster_hash_lst)) as p:
            for i_id, r_lst in p.imap_unordered(__partial_shared_simularity_recalc_entry__, iterable=parameters, chunksize=chunksize):
                for c_id, value in r_lst: 
                    item, cluster = item_lst[i_id], cluster_lst[c_id]
                    simularity_matrix[item, cluster] = min(value, 1.0)
    finally:
        global shared_hash_lst
        shared_hash_lst = None
                
    return simularity_matrix
    
def init_mp(co, arr):
    global shared_co_dct, shared_hash_lst
    shared_co_dct = co
    shared_hash_lst = arr


def calc_optimal_chunksize(size, init_chunksize) -> int:
    optimized_chunksize = min(floor(size / mp.cpu_count()), init_chunksize**2)
    optimized_chunksize = max(optimized_chunksize,1)
    return optimized_chunksize

def build_common_co_multiprocess(cluster_lst: List[Cluster], new_clusters: List[Cluster], \
    matrix: SparseDictHashMatrix[Cluster, Cluster, float], chunksize: int) ->  SparseDictHashMatrix[Cluster, Cluster, float]:


    global shared_co_dct
    if(shared_co_dct is None):
        raise Exception("shared memory co-matrix is not instantiated. call 'build_shared_memory_co_matrix' with reference to it, to instantiate it.")
    
    try:
        new_hashes = UnevenNestedList.init_ram(new_clusters)
        
        if len(cluster_lst) > 0:
            old_hashes = UnevenNestedList.init_ram(cluster_lst)
            parameters = ( (i, j, n_cls, o_cls) for i, n_cls in 
                          tqdm(enumerate(new_hashes), total=len(new_hashes), desc='adding clusters to CCO P1')\
                              for j, o_cls in enumerate(old_hashes))
            
            optimized_chunksize = calc_optimal_chunksize(len(new_hashes)*len(old_hashes), chunksize)
            #go through all new clusters and see if they are related to existing clusters
            with mp.Pool(mp.cpu_count(), initializer=init_mp, initargs=(shared_co_dct, None), maxtasksperchild=optimized_chunksize) as p:
                for id1, id2, value in p.imap_unordered(__partial_common_co__, iterable=parameters, chunksize=optimized_chunksize):
                    cluster1, cluster2 = new_clusters[id1], cluster_lst[id2]
                    matrix[cluster1, cluster2] = value
            del old_hashes
            gc.collect()
            
        #endif
        
        
        parameters = ((i, j, new_hashes[i], new_hashes[j]) for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P2')\
                            for j in range(i+1, len(new_clusters)))
        optimized_chunksize = calc_optimal_chunksize(len(new_hashes)*1.25, chunksize)
        with mp.Pool(mp.cpu_count(), initializer=init_mp, initargs=(shared_co_dct, None), maxtasksperchild=optimized_chunksize) as p:
            for id1, id2, value in p.imap_unordered(__partial_common_co__, iterable=parameters, chunksize=optimized_chunksize):
                cluster1, cluster2 = new_clusters[id1], new_clusters[id2]
                matrix[cluster1, cluster2] = value
        # print(sys.getsizeof(matrix.__internal__) )
        # print(sys.getsizeof(new_hashes))
        return matrix
    finally:
        global shared_hash_lst
        shared_hash_lst = None
    

def build_shared_memory_co_matrix(co_matrix: CoAssosiationMatrix) -> None:
    global shared_co_dct
    if shared_co_dct is not None:
        raise Exception('Shared co matrix is already instantiated. Call clear_shared_co_matrix to clear it.')    
    size = len(co_matrix)*2+1
    # keys = np.zeros(shape=size, dtype=np.int_)
    # values = np.memmap(filename=memfile, shape=size, dtype=np.double)
    
    shared_co_dct = SharedHashTable(size)
    for item_tup, value  in tqdm(co_matrix.items()):
        item1, item2 = item_tup
        index = hash((hash(item1), hash(item2)))
        shared_co_dct.insert(index, value)
        
    
    # for key, value in tqdm(tmp.items()):
    #     shared_co_dct.insert(key, value)
    shared_co_dct.test()

def clear_shared_co_matrix() -> None:
    global shared_co_dct, shared_hash_lst
    shared_co_dct = None
    shared_hash_lst = None


def __partial_shared_simularity_recalc_entry__(tuple: Tuple[object, int, List[int]]) \
    -> Tuple[object, List[Tuple[object, float]]]:
    id, item_hash, cluster_indexes= tuple 
    result = []
    
    global shared_co_dct
    
    co_matrix = shared_co_dct

    try:
        for c_index in cluster_indexes:
            value = 0.0
            cluster_hashes: List[int] = read_shared_hash_item(c_index)
            if len(cluster_hashes) == 0: continue
            for item2_hash in cluster_hashes:
                index = hash( (item_hash, int(item2_hash) ) )
                value += co_matrix.get(index, 0.0)
            value /= len(cluster_hashes) 
            result.append( (c_index, value) )
    finally:
        pass
    
    return (id, result)

def __partial_common_co__(tuple: Tuple[int, int, List[int], List[int]]) -> Tuple[int, int, float]:
    # gc.disable()
    id1, id2, c1, c2 = tuple
    global shared_co_dct
    value = __calc_common_co(c1, c2, shared_co_dct)
    return (id1, id2, value)
    


def __calc_common_co(cluster1: List[int], cluster2: List[int], co_matrix: Dict[int, float]) -> float:
    value = 0.0
    if len(cluster1) == 0 or len(cluster2) == 0: 
        return value
    for item1 in cluster1:
        for item2 in cluster2:
            entry_index = hash( (item1, item2) )
            value += co_matrix.get(entry_index, 0.0)
            
    return value / ( len(cluster1) + len(cluster2) )