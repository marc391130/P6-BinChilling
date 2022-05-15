from __future__ import annotations
import numpy as np
from tqdm import tqdm
from ClusterDomain import Cluster, PartitionSet
from typing import Tuple, Dict, List
import Assertions as Assert 
import scipy.sparse as sp  
from SparseMatrix_implementations import DoubleSparseDictHashMatrix, SortKeysByHash, SparseDictHashMatrix, SparseTupleHashMatrix
from math import sqrt
from multiprocessing import Pool, cpu_count


def cluster_simularity(cluster1: Cluster, cluster2: Cluster, total_elements: int) -> float:    
    if cluster1.SamePartitionAs(cluster2): return np.NINF
    intersection_len, len1, len2 = cluster1.intersection_len(cluster2), len(cluster1), len(cluster2)
    counter = intersection_len - ((len1 * len2) / total_elements)
    divisor = sqrt(len1 * len2 * (1 - (len1 / total_elements)) * (1 - (len2 / total_elements)))
    return min(counter / divisor, 1.0) if divisor != 0 else 0.0


class MemberSimularityMatrix(DoubleSparseDictHashMatrix[object, Cluster, float]):
    def __init__(self, item_lst: List[object], cluster_lst: List[Cluster]) -> None:
        self.all_clusters = cluster_lst
        self.all_items = item_lst
        super().__init__(default_value=0, sparse_value_predicate=lambda x: x == 0)         
    
    @staticmethod
    def IndependentBuild(clusters: List[Cluster], gamma: PartitionSet):
        cluster_lst = list(clusters)
        item_lst = list(gamma.get_all_elements().keys()) 
        
        matrix = MemberSimularityMatrix(item_lst, cluster_lst)
        # matrix = sp.csc_matrix(shape, dtype=np.float32)
        
        for cluster in cluster_lst:
            membership_map = cluster.calc_all_membersimularity(max_member_value=len(gamma))
            for item, simularity in membership_map.items():
                matrix[item, cluster] = simularity
        return matrix
    
    def set(self, __k: Tuple[object, Cluster], __v: float):
        if __v == 0.0 and self.has_tuple(__k): 
            self.pop(__k)
            return
        if(__v > 1 or __v < 0):
            raise ValueError(f"similarity value '{__v}' is outside range 0-1")
        super().set(__k, __v)
            
    def Item_sum(self, item) -> float:
        return sum(self.get_row(item).values())
    
    def item_max(self, item) -> float:
        row = self.get_row(item).values()
        if len(row) == 0: return 0.0
        return max(row)
    
    def item_mean(self, item: object) -> float:
        return self.Item_sum(item) / len(self.get_row(item))
    
    def item_argMax(self, item) -> Tuple[Cluster, float]:
        arg_max, max_value = None, -1 
        
        row = self.get_row(item)
        if row is None or len(row) == 0:
            return (None, 0)
        
        for cluster, sim in row.items():
            if sim >= max_value:
                arg_max, max_value = cluster, sim
        if arg_max is None:
            raise KeyError("Could not find index in argmax")
        return (arg_max, max_value)
    
    def cluster_sum(self, cluster: Cluster) -> float:
        return sum(self.get_column(cluster).values())
    
    def cluster_max(self, cluster: Cluster) -> float:
        return sum(self.get_column(cluster).values())
        
    def cluster_mean(self, cluster: Cluster) -> float:
        return self.cluster_sum(cluster) / len(self.get_column(cluster))
    
    def remove_cluster(self, cluster: Cluster) -> None:
        for item in self.get_column(cluster).keys():
            self.pop_entry(item, cluster)
    
    def assign_item_to(self, cluster: Cluster, item: object, update_sim: bool = True) -> None:
        for cluster2 in self.get_row(item).keys():
            self.pop_entry(item, cluster2)
        
        if update_sim:
            self.set( (item, cluster), 1)
    
    def get(self, __k: Tuple[object, Cluster]) -> float:
        if self.has_tuple(__k):
            return super().get(__k)
        return 0.0


class MemberMatrix(SparseDictHashMatrix[object, Cluster, int]):
    def __init__(self, cluster_index_map: set[Cluster], item_index_map: set[object]) -> None:
        super().__init__(default_value=0.0)
        self.items_map = item_index_map
        self.clusters = cluster_index_map
        self.common_neighbor_cache: Dict[object, List[object]] = None
        
        
    @staticmethod
    def build(cluster_lst: List[Cluster], item_lst: List[object]) -> MemberMatrix:
        item_index_map = set(item_lst)
        cluster_index_map = set(cluster_lst)
        matrix = MemberMatrix(cluster_index_map, item_index_map)
        
        for cluster in cluster_lst:
            for item, membership in cluster.calc_all_membership().items():
                matrix.set_entry(item, cluster, membership)
        
        return matrix
    
    
    # def get_all_common_items(self, totally_uncertain_items: set[object]) -> SparseDictHashMatrix[object, set[object]]:
    #     all_common_neighbors = SparseDictHashMatrix[object, set](default_value=[])

    #     for cluster in tqdm(self.clusters):
    #         intersection = totally_uncertain_items.intersection(cluster)
    #         for item1, item2 in [ (x, y) for x in intersection for y in cluster if x is not y]:
    #             if all_common_neighbors.has_entry(item1, item2) is False: all_common_neighbors[item1, item2] = set() 
    #             for item in [x for x in cluster if x is not item1 and x is not item2]:
    #                 all_common_neighbors[item1, item2].add(item)
    #                 all_common_neighbors[item2, item1].add(item)

    #     return all_common_neighbors
    
    def calculate_coassosiation_matrix(self, item_lst: set[object], gamma: PartitionSet) -> SparseDictHashMatrix[object, Tuple[float, int]]:
        all_common_neighbors = SparseDictHashMatrix[object, object, Tuple[float, int]](SortKeysByHash, default_value=(0,0))
        dct = {item: gamma.calc_all_coassosiation(item) for item in gamma.get_all_items()}

        def get(item1, item2) -> float:
            return max(\
                dct.get(item1, {}).get(item2, 0.0), \
                dct.get(item2, {}).get(item1, 0.0) \
                )

        for cluster in tqdm(self.clusters):
            intersection = item_lst.intersection(cluster)
            for item1, item2 in [ (x, y) for x in intersection for y in cluster if x is not y and y not in intersection]:
                if all_common_neighbors.has_entry(item1, item2) is False: all_common_neighbors[item1, item2] = (0, 0)
                for common_item in [x for x in cluster if x is not item1 and x is not item2]:
                    value, count = all_common_neighbors[item1, item2]
                    all_common_neighbors[item1, item2] = (value + get(item1, common_item) + get(item2, common_item), count+1)
                    
        return all_common_neighbors
    

                
                    
    
    def average_common_neighbors(self, item: object, cluster: Cluster, common_coassosiation_matrix: SparseDictHashMatrix[object, Tuple[float, int]] ) -> float:
        if len(cluster) <= 0:
            return 0
        if len(cluster) == 1 and item in cluster:
            return 0
        
        sum_value = sum([sum / (2*count) for sum, count in [common_coassosiation_matrix[item, item2] \
            for item2 in cluster if item != item2 ] if count != 0 ])
        
        return sum_value / len(cluster)
    
    def get_cluster_set(self, item_lst: List[object]) -> set[Cluster]:
        result = set()
        for item in item_lst:
            result.update(self.get_row(item).keys())
        return result
       
    
    def total_common_simularity(self, item_lst: set[object], gamma: PartitionSet)\
        -> SparseDictHashMatrix[object, object, float]:
        neighbor_simularity = self.calculate_coassosiation_matrix(item_lst, gamma)
        dct = {item: gamma.calc_all_coassosiation(item) for item in item_lst}
        result = SparseDictHashMatrix[object, float](SortKeysByHash,\
            default_value=0, sparse_value_predicate=lambda x: x == 0.0)
        
        def get(item1, item2) -> float:
            return max(\
                dct.get(item1, {}).get(item2, 0.0), \
                dct.get(item2, {}).get(item1, 0.0) \
                )
        
        for i1 in range(len(item_lst)):
            for i2 in range(i1+1, len(item_lst)):
                item1, item2 = item_lst[i1], item_lst[i2]
                result[item1, item2] = 0.5 * (get(item1, item2) + neighbor_simularity.getEntry(item1, item2) )
        return result
        

class CoAssosiationMatrix(SparseDictHashMatrix[object, object, float]):
    def __init__(self) -> None:
        super().__init__(SortKeysByHash, default_value=0.0)
        # self.matrix = matrix
        # self.index_map = index_map
        
    
    @staticmethod
    def build_old(gamma: PartitionSet) -> CoAssosiationMatrix:
        item_lst = list((gamma.get_all_elements()).keys())
        index_map = {item_lst[i]: i for i in range(len(item_lst))}
        partition_count = len(gamma)
        shape = (len(item_lst), len(item_lst))
        # matrix :np.matrix = np.full(shape=shape, fill_value=0, dtype=np.float32)
        matrix = sp.lil_matrix(shape, dtype=np.float16)
        
        for item1 in tqdm(item_lst):
            matrix_index1 = index_map[item1]
            item_coasssiation =  gamma.calc_all_coassosiation(item1)
            
            for item2, value in item_coasssiation.items():
                matrix_index2 = index_map[item2]
                    
                matrix[matrix_index1, matrix_index2] = value
                matrix[matrix_index2, matrix_index1] = value
        
        return CoAssosiationMatrix(matrix, index_map, partition_count)
    
    @staticmethod
    def build(gamma: PartitionSet, partial_lst: List[object] = None) -> CoAssosiationMatrix:
        item_lst = partial_lst if partial_lst is not None else gamma.get_all_items()
        matrix = CoAssosiationMatrix()
        
        for item1 in tqdm(item_lst):
            item_coasssiation =  gamma.calc_all_coassosiation(item1)
            for item2, value in item_coasssiation.items():
                matrix[item1, item2] = value
        
        return matrix
    
    def FindAssociatedItem(self, item: object) -> object or None:
        row = self.get_row(item)
        max_value, max_item = 0, None
        for other_item, value in row.items():
            if item is other_item: continue
            if value > max_value:
                max_value = value
                max_item = other_item
        return max_item
    
    
    def bin_mean(self, c1: Cluster, c2: Cluster) -> float:
        if len(c1) == 0 or len(c2) == 0: return 0.0        
        
        result = sum( ( 1 if x is y else self.getEntry(x, y) for x in c1 for y in c2) )
        return result / (len(c1) * len(c2))
    
    def has_association(self, c1: Cluster, c2: Cluster) -> bool:
        if len(c1) == 0 or len(c2) == 0: return False
        return any( (self.has_entry(x, y) or x is y for x in c1 for y in c2) )
    
    def cluster_mean(self, item: object, cluster: Cluster) -> float:
        if len(cluster) == 0: return 0.0
        result = sum([self.getEntry(item, x)  for x in cluster if item is not x])
        return result / len(cluster)
    

    def cluster_sim_mean(self, item: object, cluster: Cluster, simularity_matrix: MemberSimularityMatrix) -> float:
        if len(cluster) == 0: return 0.0
        result = sum([self.getEntry(item, x) * sim \
            for x, sim in simularity_matrix.get_column(item) if item is not x])
        return result / len(cluster)


class ClustserSimularityMatrix:
    def __init__(self, cluster_lst: List[Cluster], total_item_count: int, min_value: float, matrix: SparseTupleHashMatrix[Cluster, float] = None) -> None:
        self.matrix = SparseTupleHashMatrix[Cluster, float](SortKeysByHash) if matrix is None else matrix
        self.__all_clusters__ = set(cluster_lst)
        self.__non_similar_clusters__ = set(cluster_lst)
        self.total_item_count = total_item_count
        self.min_value = min_value

    @staticmethod
    def build(cluster_lst: List[Cluster], total_item_count: int, a1_min: float) -> ClustserSimularityMatrix:
        matrix = ClustserSimularityMatrix(cluster_lst, total_item_count, a1_min)
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
    def build_multithread_windows(cluster_lst: List[Cluster], total_item_count: int, a1_min: float,\
        processes: int, chunksize: int = 75) -> ClustserSimularityMatrix:
        if processes <= 1: return ClustserSimularityMatrix.build(cluster_lst, total_item_count, a1_min)
        matrix_dct = SparseTupleHashMatrix(SortKeysByHash)
        
        parameters = [(i, cluster_lst, total_item_count, a1_min) for i in range(len(cluster_lst))]
        # result : List[List[Tuple[Cluster, Cluster, float]]] = None
        with Pool(min(cpu_count(), processes)) as p:
            for chunk_result in tqdm(p.imap(partial_build_similarity_row, parameters, chunksize=chunksize), total=len(cluster_lst)):
                for i1, i2, sim in chunk_result:
                    matrix_dct[ cluster_lst[i1], cluster_lst[i2] ] = sim
        
        return ClustserSimularityMatrix(cluster_lst, total_item_count, a1_min, matrix_dct)
    
    @staticmethod
    def build_multithread(cluster_lst: List[Cluster], total_item_count: int, a1_min: float,\
        processes: int, chunksize: int = 75) -> ClustserSimularityMatrix:
        if processes <= 1: return ClustserSimularityMatrix.build(cluster_lst, total_item_count, a1_min)
        matrix_dct = SparseTupleHashMatrix(SortKeysByHash)
        try: 
            global shared_cluster_lst_memory 
            shared_cluster_lst_memory = cluster_lst
                
            parameters = [(i, total_item_count, a1_min) for i in range(len(cluster_lst))]
            # result : List[List[Tuple[Cluster, Cluster, float]]] = None
            with Pool(min(cpu_count(), processes)) as p:
                for chunk_result in tqdm(p.imap_unordered(partial_build_similarity_row_shared_mem, parameters, chunksize=chunksize), total=len(cluster_lst)):
                    for i1, i2, sim in chunk_result:
                        matrix_dct[ cluster_lst[i1], cluster_lst[i2] ] = sim
            #end pool            
            return ClustserSimularityMatrix(cluster_lst, total_item_count, a1_min, matrix_dct)
        finally:
            #del will not work here, as the variable still needs to be know, but should not hold a reference
            shared_cluster_lst_memory = None
    
    
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