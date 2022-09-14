import gc
from typing import Generator, Iterable, Iterator, List, Dict, Tuple
from Cluster_matrices import MemberSimularityMatrix
from ClusterDomain import Cluster
from multiprocessing import Pool, cpu_count
from SparseMatrix_implementations import SparseDictHashMatrix
from tqdm import tqdm
from Cluster_matrices import CoAssosiationMatrix
from sqlitedict import SqliteDict
from time import time

# (combined hash of two items) -> (co-assosiation value ( CO function ) for those two items)

def CreateDatabase(table: str) -> SqliteDict:
    return SqliteDict("co_matrix.sqlite", tablename=table,\
            flag='w',autocommit=False, journal_mode='OFF', timeout=20, outer_stack=False)

def AttachDatabase(table: str) -> SqliteDict:
    return SqliteDict("co_matrix.sqlite", tablename=table,\
            flag='c',autocommit=False, journal_mode='OFF', timeout=5, outer_stack=False)

def OpenDatabase(table: str) -> SqliteDict:
    return SqliteDict("co_matrix.sqlite", tablename=table,\
            flag='r',autocommit=False, journal_mode='OFF', outer_stack=False)

def recalculate_simularity_multiprocess(item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
    cluster_lst: List[Cluster], chunksize: int) -> MemberSimularityMatrix:
    
    #Generate the reduced cluster lst map and cluster_hash_lst
    #reduced_cluster_lst_map is a map of object, to a list of indexes for clusters
        #this is used to not have to loop over all clusters, only the relevant ones.
    reduced_cluster_lst_map: Dict[int, List[int]] = {}
    cluster_hash_lst: List[List[int]] = []
    for c_index in tqdm(range(len(cluster_lst)), desc='Prep. and reducing param.  '):
        cluster = cluster_lst[c_index]
        cluster_hash_lst.append([hash(item) for item in cluster])
        for item in cluster:
            reduced_cluster_lst_map[hash(item)] = \
                reduced_cluster_lst_map.get(hash(item), []) + [c_index]
                
    
    
    rpdb = CreateDatabase('reduced_parameters') 
    chdb = CreateDatabase('cluster_hash_lst')

    try:
        rpdb.update(reduced_cluster_lst_map.items())
        rpdb.commit(blocking=True)
        print("commited 1")
        chdb.update(enumerate(cluster_hash_lst))
        chdb.commit(blocking=True)
        print("commited 2")
        print(f"item len is: {len(item_lst)}")
        
        start = time()
        with Pool(cpu_count()) as p:
            for i_id, r_lst in tqdm(p.imap_unordered(__partial_shared_simularity_recalc_entry__, \
                iterable=((i, hash(item)) for i, item in enumerate(item_lst)), chunksize=chunksize),\
                    total=len(item_lst), desc='Recalc. similarity matrix  '):
                print(f"returned: {i_id} after {time() - start}s")
                for c_id, value in r_lst: 
                    item, cluster = item_lst[i_id], cluster_lst[c_id]
                    simularity_matrix[item, cluster] = min(value, 1.0)
        
        print("multiprocess over")
        return simularity_matrix
    finally:
        print("closing connections")
        rpdb.clear()
        rpdb.close(do_log=False)
        chdb.clear()
        chdb.close(do_log=False)

def build_common_co_multiprocess(cluster_lst: List[Cluster], new_clusters: List[Cluster], \
    matrix: SparseDictHashMatrix[Cluster, Cluster, float], chunksize: int) ->  SparseDictHashMatrix[Cluster, Cluster, float]:

    ocdb = CreateDatabase('old_clusters')
    ncdb = CreateDatabase('new_clusters')
    
    def format_clusters(cls: List[Cluster]) -> Iterable[Tuple[int, List[int]]]:
        return enumerate([hash(i) for i in c] for c in cls)
    
    try:
        ocdb.update(format_clusters(cluster_lst))
        ncdb.update(format_clusters(new_clusters))
        ocdb.commit(blocking=True)
        ncdb.commit(blocking=True)

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
        ocdb.clear()
        ocdb.close()
        ncdb.clear()
        ncdb.close()


class tqdmItemsDummy():
    def __init__(self, iterable) -> None:
        self.iterable = iterable
    def items(self):
        return ( (key, value) for key, value in tqdm(self.iterable) )

def __co_matrix_generator__(co_matrix: CoAssosiationMatrix) -> Generator:
    known_indexes = set()
    for item, index2_map  in tqdm(co_matrix.items()):
        # item, index2_map: Tuple[object, Dict[object, float]]
        for item2, value in index2_map.items():
            index = hash((hash(item), hash(item2)))
            if index in known_indexes: 
                raise Exception("non unique hash value optained from two items")
            known_indexes.add(index)
            yield (index, value)

def build_shared_memory_co_matrix(co_matrix: CoAssosiationMatrix) -> None:
    generator = __co_matrix_generator__(co_matrix)
    with CreateDatabase('co_matrix') as codb:
        codb.update(generator)
        print('now comitting')
        codb.commit(blocking=True)

def clear_shared_co_matrix() -> None:
    with AttachDatabase('co_matrix') as codb:
        codb.terminate()


def getOrDefault(db: Dict, key, default = None):
    try:
        return db[key]
    except KeyError:
        return default
    

def __partial_shared_simularity_recalc_entry__(tuple: Tuple[object, int]) \
    -> Tuple[object, List[Tuple[object, float]]]:
    id, item_hash = tuple 
    result = []
    
    co_matrix = OpenDatabase('co_matrix')
    reduced_smd = OpenDatabase('reduced_parameters')
    cluster_hash_smd = OpenDatabase('cluster_hash_lst')
    
    try:
        cluster_indexes : List[int] = getOrDefault(reduced_smd, item_hash, default=[])
        
        for c_index in cluster_indexes:
            value = 0.0
            cluster_hashes: List[int] = cluster_hash_smd[c_index]
            if len(cluster_hashes) == 0: continue
            for item2_hash in cluster_hashes:
                index = hash( (item_hash, int(item2_hash) ) )
                value += float(getOrDefault(co_matrix, index, 0.0))
            result.append( (c_index, value) )
    finally:
        co_matrix.close()
        reduced_smd.close()
        cluster_hash_smd.close()
    
    return (id, result)
        
    


def __partial_shared_oldnew_common_co_entry(tuple: Tuple[object, List[int]]):
    id, cluster1 = tuple
    if len(cluster1) == 0: return (id, [])
    result = []
    
    co_matrix = OpenDatabase('co_matrix')
    old_hash_smd = OpenDatabase('old_clusters')
    
    try:
        for index, cluster2 in old_hash_smd.items():
            if len(cluster2) == 0: continue 
            value = __calc_common_co(cluster1, cluster2, co_matrix)
            result.append( (int(index), value) )
            
        return (id, result)
    finally:
        co_matrix.close()
        old_hash_smd.close()


def __partial_shared_newnew_common_co_entry(tuple: Tuple[object, int]):
    id, index = tuple
    result = []
    
    co_matrix = OpenDatabase('co_matrix')
    new_hash_smd = OpenDatabase('new_clusters')
    
    try:
        hash_cluster: List[int] = new_hash_smd[index]
        if len(hash_cluster) == 0: return (id, [])
        for index in range(index + 1, len(new_hash_smd)):
            cluster2 = new_hash_smd[index]
            if len(cluster2) == 0: continue
            value = __calc_common_co(hash_cluster, cluster2, co_matrix)
            result.append( (int(index), value) )
            
        return (id, result)
    finally:
        co_matrix.close()
        new_hash_smd.close()
        

def __calc_common_co(cluster1: List[int], cluster2: List[int], co_matrix: Dict[int, float]) -> float:
    value = 0.0
    if len(cluster1) == 0 or len(cluster2) == 0: 
        return value
    for item1 in cluster1:
        for item2 in cluster2:
            entry_index = hash( (item1, item2) )
            value += float(getOrDefault(co_matrix, entry_index, 0.0))
            
    return value / ( len(cluster1) + len(cluster2) )