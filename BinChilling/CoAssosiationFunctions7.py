from multiprocessing.managers import DictProxy, SyncManager
from sys import getsizeof
from typing import List, Dict, Tuple
from Cluster_matrices import MemberSimularityMatrix
from ClusterDomain import Cluster
import multiprocessing as mp
from SparseMatrix_implementations import SparseDictHashMatrix
from tqdm import tqdm
from Cluster_matrices import CoAssosiationMatrix
import numpy as np

# (combined hash of two items) -> (co-assosiation value ( CO function ) for those two items)
shared_hash_lst: List[List[int]] = None

def read_shared_hash_item(index: int) -> List[int]:
    global shared_hash_lst
    return shared_hash_lst[index]

def recalculate_simularity_multiprocess(item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
    cluster_lst: List[Cluster], chunksize: int) -> MemberSimularityMatrix:
    
    #Generate the reduced cluster lst map and cluster_hash_lst
    #reduced_cluster_lst_map is a map of object, to a list of indexes for clusters
        #this is used to not have to loop over all clusters, only the relevant ones.
    reduced_cluster_lst_map: Dict[object, List[int]] = {}
    cluster_hash_lst: List[str] = np.empty(len(cluster_lst), dtype=object)
    for c_index in tqdm(range(len(cluster_lst)), desc='Prep. and reducing param.  '):
        cluster = cluster_lst[c_index]
        cluster_hash_lst[c_index] = np.array([hash(item) for item in cluster], dtype=np.int_)
        for item in cluster:
            reduced_cluster_lst_map[item] = \
                reduced_cluster_lst_map.get(item, []) + [c_index]
        
    try:
        global manager_proxy
        manager, proxy = manager_proxy
        manager_proxy = None #dont share from fork
    
        parameters = ((i, hash(item), reduced_cluster_lst_map.get(item, []), proxy) for i, item in tqdm(enumerate(item_lst), desc='Recalc. similarity matrix  ', total=len(item_lst)))

        
        with mp.Pool(mp.cpu_count(), initializer=init_mp, initargs=(cluster_hash_lst, )) as p:
            for i_id, r_lst in p.imap_unordered(__partial_shared_simularity_recalc_entry__, iterable=parameters, chunksize=chunksize):
                for c_id, value in r_lst: 
                    item, cluster = item_lst[i_id], cluster_lst[c_id]
                    simularity_matrix[item, cluster] = min(value, 1.0)
    finally:
        manager_proxy = (manager, proxy)
                
    return simularity_matrix
    
def init_mp(arr):
    global shared_hash_lst
    shared_hash_lst = arr

def build_common_co_multiprocess(cluster_lst: List[Cluster], new_clusters: List[Cluster], \
    matrix: SparseDictHashMatrix[Cluster, Cluster, float], chunksize: int) ->  SparseDictHashMatrix[Cluster, Cluster, float]:

    def format_clusters(cls: List[Cluster]) -> List[List[int]]:
        return np.array([format_cls(c) for c in cls], dtype=object)
    
    def format_cls(c: Cluster) -> List[int]:
        return np.array([hash(i) for i in c], dtype=np.int_)
    
    new_hashes = format_clusters(new_clusters)
    
    try:
        global manager_proxy
        manager, proxy = manager_proxy
        manager_proxy = None #dont share from fork
        
        if len(cluster_lst) > 0:
            parameters = ( (i, enc_cls) for i, enc_cls in tqdm(enumerate(new_hashes), total=len(new_hashes), desc='adding clusters to CCO P1') )
            old_hashes = format_clusters(cluster_lst)
            #go through all new clusters and see if they are related to existing clusters
            with mp.Pool(mp.cpu_count(), initializer=init_mp, initargs=(old_hashes, )) as p:
                for id1, r_val in p.imap_unordered(__partial_shared_oldnew_common_co_entry, iterable=parameters, chunksize=chunksize):
                    for id2, value in r_val:
                        cluster1, cluster2 = new_clusters[id1], cluster_lst[id2]
                        matrix[cluster1, cluster2] = value
        #endif
        
        parameters = ( (i, i) for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P2') )
        with mp.Pool(mp.cpu_count(), initializer=init_mp, initargs=(new_hashes, )) as p:
            for id1, r_val in p.imap_unordered(__partial_shared_newnew_common_co_entry, iterable=parameters, chunksize=chunksize):
                for id2, value in r_val:  
                    cluster1, cluster2 = new_clusters[id1], new_clusters[id2]
                    matrix[cluster1, cluster2] = value
        return matrix
    finally:
        manager_proxy = (manager, proxy)


manager_proxy: Tuple[SyncManager, DictProxy] = None
def build_shared_memory_co_matrix(co_matrix: CoAssosiationMatrix) -> None:
    global manager_proxy
    dic = {}
    for item, index2_map  in tqdm(co_matrix.items()):
        # item, index2_map: Tuple[object, Dict[object, float]]
        for item2, value in index2_map.items():
            index = hash((hash(item), hash(item2)))
            if index in dic: 
                raise Exception("non unique hash value optained from two items")
            dic[index] = value
    print(getsizeof(dic))
    manager = mp.Manager()
    global_dct = manager.dict()
    global_dct.update(dic)
    manager_proxy = (manager, global_dct)
    

def clear_shared_co_matrix() -> None:
    global manager_proxy
    manager, global_dct = manager_proxy
    manager.shutdown()


def __partial_shared_simularity_recalc_entry__(tuple: Tuple[object, int, List[int], Dict[int, float]]) \
    -> Tuple[object, List[Tuple[object, float]]]:
    id, item_hash, cluster_indexes, co_matrix= tuple 
    result = []
    
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
        
    


def __partial_shared_oldnew_common_co_entry(tuple: Tuple[object, List[int]]):
    id, cluster1 = tuple
    if len(cluster1) == 0: return (id, [])
    result = []
    
    global shared_co_dct, shared_hash_lst
    co_matrix = shared_co_dct
    
    try:
        for index in range(len(shared_hash_lst)):
            cluster2 = read_shared_hash_item(index)
            if len(cluster2) == 0: continue 
            value = __calc_common_co(cluster1, cluster2, co_matrix)
            result.append( (int(index), value) )
            
        return (id, result)
    finally:
        pass


def __partial_shared_newnew_common_co_entry(tuple: Tuple[object, int]):
    id, index = tuple
    result = []
    
    global shared_co_dct, shared_hash_lst
    co_matrix = shared_co_dct
    
    new_hash_smd = shared_hash_lst
    
    try:
        hash_cluster: List[int] = read_shared_hash_item(index)
        if len(hash_cluster) == 0: return (id, [])
        for index in range(index + 1, len(new_hash_smd)):
            cluster2 = read_shared_hash_item(index)
            if len(cluster2) == 0: continue
            value = __calc_common_co(hash_cluster, cluster2, co_matrix)
            result.append( (int(index), value) )
            
        return (id, result)
    finally:
        pass
        

def __calc_common_co(cluster1: List[int], cluster2: List[int], co_matrix: Dict[int, float]) -> float:
    value = 0.0
    if len(cluster1) == 0 or len(cluster2) == 0: 
        return value
    for item1 in cluster1:
        for item2 in cluster2:
            entry_index = hash( (item1, item2) )
            value += co_matrix.get(entry_index, 0.0)
            
    return value / ( len(cluster1) + len(cluster2) )