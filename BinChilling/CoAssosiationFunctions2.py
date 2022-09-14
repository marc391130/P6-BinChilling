import gc
from typing import Iterator, List, Dict, Tuple
from Cluster_matrices import MemberSimularityMatrix
from ClusterDomain import Cluster
from multiprocessing import Pool, cpu_count, Array
from SparseMatrix_implementations import SparseDictHashMatrix
from tqdm import tqdm
from Cluster_matrices import CoAssosiationMatrix
import shared_memory_dict as smd

# (combined hash of two items) -> (co-assosiation value ( CO function ) for those two items)

shared_len : Tuple[int, int] = None

def recalculate_simularity_multiprocess(item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
    cluster_lst: List[Cluster], chunksize: int) -> MemberSimularityMatrix:
    
    global shared_len
    #Generate the reduced cluster lst map and cluster_hash_lst
    #reduced_cluster_lst_map is a map of object, to a list of indexes for clusters
        #this is used to not have to loop over all clusters, only the relevant ones.
    reduced_cluster_lst_map: Dict[object, List[int]] = {}
    cluster_hash_lst: List[List[int]] = []
    for c_index in tqdm(range(len(cluster_lst)), desc='Prep. and reducing param.  '):
        cluster = cluster_lst[c_index]
        cluster_hash_lst.append([hash(item) for item in cluster])
        for item in cluster:
            reduced_cluster_lst_map[item] = \
                reduced_cluster_lst_map.get(item, []) + [c_index]
                
    reduced_smd = smd.SharedMemoryDict('reduced_parameters', shared_memory_size(len(reduced_cluster_lst_map)))
    cluster_hash_smd = smd.SharedMemoryDict('cluster_hash_lst', shared_memory_size(len(cluster_hash_lst)))
    shared_len = ( len(reduced_cluster_lst_map), len(cluster_hash_lst) )
    try:
        
        __update_shared_memory_fast(reduced_smd, reduced_cluster_lst_map.items() )
        __update_shared_memory_fast(cluster_hash_smd, enumerate(cluster_hash_lst))
            
        parameters = ((i, hash(item)) for i, item in tqdm(enumerate(item_lst), desc='Recalc. similarity matrix  '))
        
        with Pool(cpu_count()) as p:
            for i_id, r_lst in p.imap_unordered(__partial_shared_simularity_recalc_entry__, iterable=parameters, chunksize=chunksize):
                for c_id, value in r_lst: 
                    item, cluster = item_lst[i_id], cluster_lst[c_id]
                    simularity_matrix[item, cluster] = min(value, 1.0)
                    gc.collect()
        
        return simularity_matrix
    finally:
        shared_len = None
        reduced_smd.shm.close()
        reduced_smd.shm.unlink()
        cluster_hash_smd.shm.close()
        cluster_hash_smd.shm.unlink()
        gc.collect()

def build_common_co_multiprocess(cluster_lst: List[Cluster], new_clusters: List[Cluster], \
    matrix: SparseDictHashMatrix[Cluster, Cluster, float], chunksize: int) ->  SparseDictHashMatrix[Cluster, Cluster, float]:

    global shared_len

    old_hash_smd = smd.SharedMemoryDict('old_clusters', shared_memory_size(len(cluster_lst)))    
    new_hash_smd = smd.SharedMemoryDict('new_clusters', shared_memory_size(len(new_clusters)))
    shared_len = ( len(cluster_lst), len(new_clusters) )
    
    try:
        for i, o_cluster in enumerate(cluster_lst):
            old_hash_smd[i] = [hash(item) for item in o_cluster]
        for i, n_cluster in enumerate(new_clusters):
            new_hash_smd[i] = [hash(item) for item in n_cluster]    

        if len(cluster_lst) > 0:
            parameters = ( i for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P1') )
        
            #go through all new clusters and see if they are related to existing clusters
            with Pool(cpu_count()) as p:
                for id1, r_val in p.imap_unordered(__partial_shared_oldnew_common_co_entry, iterable=parameters, chunksize=chunksize):
                    for id2, value in r_val:
                        cluster1, cluster2 = new_clusters[id1], cluster_lst[id2]
                        matrix[cluster1, cluster2] = value
        #endif
        
        parameters = ( (i, i) for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P2') )
        with Pool(cpu_count()) as p:
            for id1, r_val in p.imap_unordered(__partial_shared_newnew_common_co_entry, iterable=parameters, chunksize=chunksize):
                for id2, value in r_val:  
                    cluster1, cluster2 = new_clusters[id1], cluster_lst[id2]
                    matrix[cluster1, cluster2] = value
        return matrix
    finally:
        shared_len = None
        old_hash_smd.shm.close()
        old_hash_smd.shm.unlink()
        new_hash_smd.shm.close()
        new_hash_smd.shm.unlink()
        gc.collect()


smd_len = 0
def build_shared_memory_co_matrix(co_matrix: CoAssosiationMatrix) -> None:
    global smd_len
    smd_len = co_matrix.count_entries()
    shared_co_dct = smd.SharedMemoryDict('co_matrix', shared_memory_size(smd_len))
    local_dct = {}
    try:
        for item, index2_map  in tqdm(co_matrix.items()):
            # item, index2_map: Tuple[object, Dict[object, float]]
            for item2, value in index2_map.items():
                index = hash((hash(item), hash(item2)))
                if index in local_dct: 
                    raise Exception("non unique hash value optained from two items")
                local_dct[index] = value
        #end for
        __update_shared_memory_fast(shared_co_dct, local_dct.items())
    finally:
        shared_co_dct.shm.close()

def shared_memory_size(length: int) -> int:
    return length*1024*1024

def __update_shared_memory_fast(smd_dict: smd.SharedMemoryDict, values: Iterator) -> smd.SharedMemoryDict:
    with smd_dict._modify_db() as db:
        for key, value in values:
            db[key] = value
        


def clear_shared_co_matrix() -> None:
    global smd_len
    shared_co_dct = smd.SharedMemoryDict('co_matrix', shared_memory_size(smd_len))
    smd_len = 0
    shared_co_dct.shm.close()
    shared_co_dct.shm.unlink()
    gc.collect()


def __partial_shared_simularity_recalc_entry__(tuple: Tuple[object, int]) \
    -> Tuple[object, List[Tuple[object, float]]]:
    id, item_hash = tuple 
    result = []
    global smd_len, shared_len
    
    co_matrix = smd.SharedMemoryDict('co_matrix', shared_memory_size(smd_len))
    reduced_smd = smd.SharedMemoryDict('reduced_parameters', shared_memory_size(shared_len[0]))
    cluster_hash_smd = smd.SharedMemoryDict('cluster_hash_lst', shared_memory_size(shared_len[1]))
    
    try:
        cluster_indexes : List[int] = reduced_smd[item_hash]
        
        for c_index in cluster_indexes:
            value = 0.0
            cluster_hashes: List[int] = cluster_hash_smd[c_index]
            if len(cluster_hashes) == 0: continue
            for item2_hash in cluster_hashes:
                index = hash( (item_hash, int(item2_hash) ) )
                value += co_matrix.get(index, 0.0)
            result.append( (c_index, value) )
    finally:
        co_matrix.shm.close()
        reduced_smd.shm.close()
        cluster_hash_smd.shm.close()
    
    return (id, result)
        
    


def __partial_shared_oldnew_common_co_entry(tuple: Tuple[object, List[int]]):
    id, cluster1 = tuple
    if len(cluster1) == 0: return (id, [])
    
    global smd_len, shared_len
    old_len, _ = shared_len
    result = []
    
    co_matrix = smd.SharedMemoryDict('co_matrix', shared_memory_size(smd_len))
    old_hash_smd = smd.SharedMemoryDict('old_clusters', shared_memory_size(old_len))
    
    try:
        for index, cluster2 in old_hash_smd.items():
            if len(cluster2) == 0: continue 
            value = __calc_common_co(cluster1, cluster2, co_matrix)
            result.append( (int(index), value) )
            
        return (id, result)
    finally:
        co_matrix.shm.close()
        old_hash_smd.shm.close()


def __partial_shared_newnew_common_co_entry(tuple: Tuple[object, int]):
    id, index = tuple
    global smd_len, shared_len
    _, new_len = shared_len
    result = []
    
    co_matrix = smd.SharedMemoryDict('co_matrix', shared_memory_size(smd_len))
    new_hash_smd = smd.SharedMemoryDict('new_clusters', shared_memory_size(new_len))
    
    try:
        hash_cluster: List[int] = new_hash_smd[index]
        if len(hash_cluster) == 0: return (id, [])
        for index in range(index + 1, new_len):
            cluster2 = new_hash_smd[index]
            if len(cluster2) == 0: continue
            value = __calc_common_co(hash_cluster, cluster2, co_matrix)
            result.append( (int(index), value) )
            
        return (id, result)
    finally:
        co_matrix.shm.close()
        new_hash_smd.shm.close()

def __calc_common_co(cluster1: List[int], cluster2: List[int], co_matrix: Dict[int, float]) -> float:
    value = 0.0
    if len(cluster1) == 0 or len(cluster2) == 0: 
        return value
    for item1 in cluster1:
        for item2 in cluster2:
            entry_index = hash( (item1, item2) )
            value += float(co_matrix.get(entry_index, 0.0))
            
    return value / ( len(cluster1) + len(cluster2) )