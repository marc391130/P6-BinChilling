from typing import List, Dict, Tuple
from Cluster_matrices import MemberSimularityMatrix
from ClusterDomain import Cluster
from multiprocessing import Pool, cpu_count
from SparseMatrix_implementations import SparseDictHashMatrix
from tqdm import tqdm
from Cluster_matrices import CoAssosiationMatrix
import itertools

# (combined hash of two items) -> (co-assosiation value ( CO function ) for those two items)
shared_co_dct: Dict[int, float] = None


def recalculate_simularity_multiprocess(item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
    cluster_lst: List[Cluster], chunksize: int) -> MemberSimularityMatrix:
    
    global shared_co_dct
    if(shared_co_dct is None):
        raise Exception("shared memory co-matrix is not instantiated. call 'build_shared_memory_co_matrix' with reference to it, to instantiate it.")
    
    
    #Generate the reduced cluster lst map and cluster_hash_lst
    #reduced_cluster_lst_map is a map of object, to a list of indexes for clusters
        #this is used to not have to loop over all clusters, only the relevant ones.
        #
    reduced_cluster_lst_map: Dict[object, List[int]] = {}
    cluster_hash_lst: List[List[int]] = []
    for c_index in tqdm(range(len(cluster_lst)), desc='Prep. and reducing param.  '):
        cluster = cluster_lst[c_index]
        cluster_hash_lst.append([hash(item) for item in cluster])
        for item in cluster:
            reduced_cluster_lst_map[item] = \
                reduced_cluster_lst_map.get(item, []) + [c_index]
    
    parameters = (((x, c), hash(item_lst[x]), cluster_hash_lst[c]) for x in tqdm(range(len(item_lst)), desc='Recalc. similarity matrix  ')\
        for c in reduced_cluster_lst_map.get(item_lst[x], []))
    
    with Pool(cpu_count()) as p:
        for id, value in p.imap_unordered(__partial_simularity_recalc_entry__, iterable=parameters, chunksize=chunksize):
            id1, id2 = id
            item, cluster = item_lst[id1], cluster_lst[id2]
            simularity_matrix[item, cluster] = min(value, 1.0)
            
    return simularity_matrix
    

def build_common_co_multiprocess(cluster_lst: List[Cluster], new_clusters: List[Cluster], \
    matrix: SparseDictHashMatrix[Cluster, Cluster, float], chunksize: int) ->  SparseDictHashMatrix[Cluster, Cluster, float]:
    
    global shared_co_dct
    if(shared_co_dct is None):
        raise Exception("shared memory co-matrix is not instantiated. call 'build_shared_memory_co_matrix' with reference to it, to instantiate it.")
    
    new_hash_lst = [[hash(item) for item in cluster] for cluster in new_clusters]
    
    
    if len(cluster_lst) > 0:
        old_hash_lst = [[hash(item) for item in cluster] for cluster in cluster_lst]
        parameters = ( ((i, j), new_hash_lst[i], old_hash_lst[j]) \
            for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P1') \
                for j in range(len(cluster_lst)))
        
        #go through all new clusters and see if they are related to existing clusters
        with Pool(cpu_count()) as p:
            for id, value in p.imap_unordered(__partial_common_co_entry__, iterable=parameters, chunksize=chunksize):
                id1, id2 = id
                cluster1, cluster2 = new_clusters[id1], cluster_lst[id2]
                matrix[cluster1, cluster2] = value
        del old_hash_lst, parameters
    #end if
    
    #new go through all the new clusters to see if they are related
    parameters2 = ( ((i, j), new_hash_lst[i], new_hash_lst[j]) \
        for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P2') \
            for j in range(i+1, len(new_clusters)))
    
    with Pool(cpu_count()) as p:
        for id, value in p.imap_unordered(__partial_common_co_entry__, iterable=parameters2, chunksize=chunksize):
            id1, id2 = id
            cluster1, cluster2 = new_clusters[id1], new_clusters[id2]
            matrix[cluster1, cluster2] = value      
    
    return matrix


def build_shared_memory_co_matrix(co_matrix: CoAssosiationMatrix) -> None:
    global shared_co_dct
    if shared_co_dct is not None:
        raise Exception('Shared co matrix is already instantiated. Call clear_shared_co_matrix to clear it.')    
    shared_co_dct = {}
    for item, index2_map  in tqdm(co_matrix.items()):
        # item, index2_map: Tuple[object, Dict[object, float]]
        for item2, value in index2_map.items():
            index = hash((hash(item), hash(item2)))
            if index in shared_co_dct: 
                raise Exception("non unique hash value optained from two items")
            shared_co_dct[index] = value

def clear_shared_co_matrix() -> None:
    global shared_co_dct
    shared_co_dct = None


def __partial_simularity_recalc_entry__(tuple: Tuple[object, int, List[int]]) \
    -> Tuple[object, float]:
    id, item_id, cluster = tuple 
    global shared_co_dct
    value = 0.0
    
    if len(cluster) == 0:
        return (id, value)
    
    value = sum((shared_co_dct.get(hash((item_id, o_item)), 0.0) for o_item in cluster))
    
    value = value / len(cluster)
    return (id, value)

def __partial_common_co_entry__(tuple: Tuple[object, List[int], List[int]]) \
    -> Tuple[object, float]:
    id, cluster1, cluster2 = tuple
    global shared_co_dct
    value = 0.0
    
    for item1 in cluster1:
        for item2 in cluster2:
            index = hash((item1, item2))
            if index in shared_co_dct:
                value += shared_co_dct[index]
    
    value = value / (len(cluster1) + len(cluster2))
    return (id, value)