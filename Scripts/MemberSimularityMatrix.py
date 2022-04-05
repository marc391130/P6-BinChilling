from __future__ import annotations
from operator import index
import numpy as np
from sklearn import cluster
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from typing import Iterable, Tuple, Dict, List
import Assertions as Assert 

def average(ite :Iterable) -> float:
    return sum(ite) / len(ite)


class AceMatrix:
    def __init__(self) -> None:
        pass
    

class MemberSimularityMatrix:
    def __init__(self, matrix: np.matrix, cluster_index_map: Dict[Cluster, int], item_index_map: Dict[Cluster, int]) -> None:
        self.matrix = matrix
        self.cluster_index_map = cluster_index_map
        self.item_index_map = item_index_map
    
    @staticmethod
    def Build(memberMatrix: MemberMatrix) -> MemberSimularityMatrix:
        cluster_index_map = dict(memberMatrix.cluster_index_map)
        item_index_map = dict(memberMatrix.item_index_map)
        shape = (len(item_index_map), len(cluster_index_map)) 
        matrix = np.full(shape=shape, fill_value=-1, dtype=np.float32)
        
        def row_divisor(item_index) -> float:
            # return np.max(memberMatrix.matrix[item_index])
            return np.sum(memberMatrix.matrix[item_index])
        
        for item, e_index in item_index_map.items():
            row_value = row_divisor(e_index)
            for cluster, c_index in cluster_index_map.items():
                membership_value = memberMatrix[item, cluster] 
                simularityValue = membership_value / row_value if row_value > 0 else 0
                matrix[e_index, c_index] = simularityValue
        np.savetxt('simularity_matrix.txt', matrix, fmt='%f', delimiter=', ', newline='\n\n')
        return MemberSimularityMatrix(matrix, cluster_index_map, item_index_map) 
    
    
    @staticmethod
    def RefineTo(simularity_matrix: MemberSimularityMatrix, keep_clusters: List[Cluster]) -> MemberSimularityMatrix:
        Assert.assert_unique(keep_clusters)
        for keep_cluster in keep_clusters:
            Assert.assert_key_exists(keep_cluster, simularity_matrix.cluster_index_map)
        keep_cluster_index_map = {keep_clusters[i]: i for i in range(len(keep_clusters)) }           
        item_index_map = dict(simularity_matrix.item_index_map)
        shape = (len(item_index_map), len(keep_cluster_index_map))
        
        matrix = np.full(shape=shape, fill_value=0, dtype=np.float32)
        
        for item, i_index in item_index_map.items():
            for cluster, c_index in keep_cluster_index_map.items():
                matrix[i_index, c_index] = simularity_matrix[item, cluster]
        np.savetxt('simularity_matrix_refined.txt', matrix, fmt='%f', delimiter=', ', newline='\n\n')
        return MemberSimularityMatrix(matrix, keep_cluster_index_map, item_index_map)
    
    
    def shape(self) -> Tuple[int, int]:
        return (len(self.cluster_index_map), len(self.item_index_map))
    
    def __getitem__(self, tuple: Tuple[object, Cluster]) -> float:
        item, cluster = tuple
        return self.GetEntry(cluster, item)
    
    def __setitem__(self, __k: Tuple[object, Cluster], __v: float):
        if(__v > 1 or __v < 0):
            raise ValueError("similarity value cannot be outside range 0-1")
        item, cluster = __k
        i_index, c_index = self.item_index_map[item], self.cluster_index_map[cluster]
        self.matrix[i_index, c_index] = __v
    
    
    def Cluster_Mean(self, cluster: Cluster) -> float:
        Assert.assert_key_exists(cluster, self.cluster_index_map)
        if len(cluster) == 0:
            return 0.0
        return self.Cluster_Sum(cluster) / len(cluster)
    
    def Cluster_Max(self, cluster: Cluster) -> float:
        c_index, max_value = self.cluster_index_map[cluster], 0
        for item in cluster:
            max_value = max(self.matrix[self.item_index_map[item], c_index], max_value)
        return max_value
        
    def Cluster_Sum(self, cluster: Cluster) -> float:
        c_index, sum_value = self.cluster_index_map[cluster], 0
        for item in cluster:
            i_index = self.item_index_map[item]
            sum_value += self.matrix[i_index, c_index]
        return sum_value
    
    def Item_sum(self, item) -> float:
        e_index, sum_value = self.item_index_map[item], 0
        for cluster, c_index in self.cluster_index_map.items():
            sum_value += self.matrix[e_index, c_index]
        return sum_value
    
    def item_max(self, item) -> float:
        e_index = self.item_index_map[item]
        return np.max(self.matrix[e_index])
    
    def item_argMax(self, item) -> Cluster:
        e_index = self.item_index_map[item]
        index = np.argmax(self.matrix[e_index])
        
        for cluster, c_index in self.cluster_index_map.items():
            if c_index == index:
                return cluster
        raise KeyError("Could not find index in argmax")
    
    # def Remove_Clusters(self, clusters_to_remove: List[Cluster]) -> None:
    #     new_cluster_map = dict(self.cluster_index_map)
    #     for cluster in clusters_to_remove:
    #         Assert.assert_key_exists(cluster, new_cluster_map)
    #         new_cluster_map.pop(cluster)
    #     return
    
    def assign_item_to(self, cluster: Cluster, item: object) -> None:
        i_index = self.item_index_map[item]
        keep_index = self.cluster_index_map[cluster]
        for c_index in self.cluster_index_map.values():
            if keep_index == c_index:
                continue
            self.matrix[i_index, c_index] = 0
    
    
    def GetEntry(self, cluster: Cluster, item: object) -> float:
        Assert.assert_key_exists(cluster, self.cluster_index_map)
        Assert.assert_key_exists(item, self.item_index_map)
        c_index, e_index = self.cluster_index_map[cluster], self.item_index_map[item]
        return self.matrix[e_index, c_index]

class MemberMatrix:
    def __init__(self, matrix: np.matrix, cluster_index_map: Dict[Cluster, int], item_index_map: Dict[object, int]) -> None:
        self.item_index_map = item_index_map
        self.cluster_index_map = cluster_index_map
        self.matrix = matrix
        
        
    @staticmethod
    def build(cluster_lst: List[Cluster], item_lst: List[object]) -> MemberMatrix:
        item_index_map = {item_lst[i]: i for i in range(len(item_lst))}
        cluster_index_map = {cluster_lst[i]: i for i in range(len(cluster_lst))}
        shape = (len(item_index_map), len(cluster_index_map))
        matrix: np.matrix = np.zeros(shape=shape, dtype=np.int)
        
        for cluster, c_index in cluster_index_map.items():
            value_map = cluster.calc_all_membership()
            for item, value in value_map.items():
                Assert.assert_not_equal(value, 0)
                e_index = item_index_map[item]
                matrix[e_index, c_index] = value
        
        # for row in matrix:
        #     print(row.sum())
        #     Assert.assert_equal(row.sum(), 5)
        np.savetxt("member_matrix_sum_0.txt", matrix.sum(axis=0), fmt='%d', delimiter=',', newline='\n\n')
        np.savetxt("member_matrix_sum_1.txt", matrix.sum(axis=1), fmt='%d', delimiter=',', newline='\n\n')
        np.savetxt("member_matrix.txt", matrix, fmt='%d', delimiter=',', newline='\n\n')
        return MemberMatrix(matrix, cluster_index_map, item_index_map)
    
    @staticmethod
    def buildButGay(cluster_lst: List[Cluster], item_lst: List[object]) -> MemberMatrix:
        item_index_map = {item_lst[i]: i for i in range(len(item_lst))}
        cluster_index_map = {cluster_lst[i]: i for i in range(len(cluster_lst))}
        shape = (len(item_index_map), len(cluster_index_map))
        matrix: np.matrix = np.zeros(shape=shape, dtype=np.int)
        
        print("Calculating gay matrix")
        for item, i_index in tqdm(item_index_map.items()):
            for cluster, c_index in cluster_index_map.items():
                value = cluster.calc_membership(item)
                matrix[i_index, c_index] = value
        
        for row in matrix:
            print(row.sum())
            Assert.assert_equal(row.sum(), 5)
        
        gay = MemberMatrix(matrix, cluster_index_map, item_index_map)
        old = MemberMatrix.build(cluster_lst, item_lst)
        
        print("Determining gayness")
        for item, i_index in tqdm(item_index_map.items()):
            for cluster, c_index in cluster_index_map.items():
                Assert.assert_equal(gay.matrix[i_index, c_index], old.matrix[i_index, c_index])
        
        print("Not gay")
        return MemberMatrix(matrix, cluster_index_map, item_index_map)
    
    
    def __getitem__(self, tuple: Tuple[object, Cluster]) -> float:
        item, cluster = tuple
        return self.getEntry(cluster, item)
     
    def getEntry(self, cluster: Cluster, item: object) -> int:
        Assert.assert_key_exists(cluster, self.cluster_index_map)
        Assert.assert_key_exists(item, self.item_index_map)
        c_index = self.cluster_index_map[cluster]
        e_index = self.item_index_map[item]
        return self.matrix[e_index, c_index]
        
    def __len__(self) -> int:
        return len(self.cluster_index_map)
    
    
    def get_cluster_row(self, cluster: Cluster) -> Dict[object, int]:
        Assert.assert_key_exists(cluster, self.cluster_index_map)
        c_index = self.cluster_index_map[cluster]
        return {item: self.matrix[e_index, c_index] for item, e_index in self.item_index_map.items()}
    
    def add_cluster(self, cluster: Cluster) -> None:
        item_map = cluster.calc_all_membership()
        column = [[item_map[item] if item in item_map else 0] for item, i in self.item_index_map.items()]
        self.cluster_index_map[cluster] = len(self.cluster_index_map)
        
        self.matrix = np.append(self.matrix, column, axis=1)
    
    def BuildSimularityMatrix(self, clusters: List[Cluster]) -> MemberSimularityMatrix:
        for cluster in clusters: # make sure every cluster exists inside the cluster_index_map
            Assert.assert_key_exists(cluster, self.cluster_index_map)
        
        for cluster, c_index in self.cluster_index_map.items():
            item_map = cluster.calc_all_membership()
            for item, i_index in self.item_index_map.items():
                if item in item_map:
                    Assert.assert_equal(self.matrix[i_index, c_index], item_map[item])
                else:
                    Assert.assert_equal(self.matrix[i_index, c_index], 0)
        return MemberSimularityMatrix.Build(self)       
    
    
    #returns a set of items
    def get_common_items(self, item1: object, item2: object) -> set:
        Assert.assert_index_exists(item1, self.item_index_map)
        Assert.assert_index_exists(item2, self.item_index_map)
        Assert.assert_not_equal(item1, item2)
        i1, i2 = self.item_index_map[item1], self.item_index_map[item2]
        Assert.assert_not_equal(i1, i2)
        
        row1 = self.matrix[i1]
        row2 = self.matrix[i2]
        Assert.assert_equal(len(row1), len(row2))
        Assert.assert_not_equal(len(row1), 0)
        common_index = set([index for index in range(len(row1)) if row1[index] > 0 and row2[index] > 0])
        common_neighbors = set()
        for cluster, c_index in self.cluster_index_map.items():
            if c_index not in common_index:
                continue
            for item in cluster:
                if item not in common_neighbors:
                    common_neighbors.add(item)
        return common_neighbors
    
    def common_neighbors(self, coassociation_matrix: CoAssosiationMatrix, item1: object, item2: object) -> float:
        common_items = self.get_common_items(item1, item2)
        sum_value = sum([coassociation_matrix[item1, a_item] + coassociation_matrix[a_item, item2] for a_item in common_items])
        div = 2* len(common_items)
        return sum_value / div if div != 0 else 0

    def average_common_neighbors(self, coassociation_matrix: CoAssosiationMatrix, item: object, cluster: Cluster):
        if len(cluster) <= 0:
            return 0
        if len(cluster) == 1 and item in cluster:
            return 0
        values = [self.common_neighbors(coassociation_matrix, item, item2) for item2 in cluster if item != item2 ]
        return sum(values) / len(values)
        
        
    def buildCoAssosiationMatrix(self, gamma: PartitionSet) -> CoAssosiationMatrix:
        return CoAssosiationMatrix.build(gamma, list(self.cluster_index_map.keys()))
    
    def shape(self) -> Tuple[int, int]:
        return (len(self.item_index_map), len(self.cluster_index_map))

class CoAssosiationMatrix:
    def __init__(self, matrix: np.matrix, index_map: Dict[object, int], partition_count: int) -> None:
        self.matrix = matrix
        self.index_map = index_map
        self.__parition_count__ = partition_count
        
    
    @staticmethod
    def build(gamma: PartitionSet, cluster_lst: List[Cluster]) -> CoAssosiationMatrix:
        Assert.assert_unique(cluster_lst)
        item_lst = list((gamma.get_all_elements()).keys())
        index_map = {item_lst[i]: i for i in range(len(item_lst))}
        partition_count = len(gamma)
        
        matrix = np.full(shape=(len(item_lst), len(item_lst)), fill_value=0, dtype=np.float32)
        
        print("Building coassociation matrix...")
        for item_index1 in tqdm(range(len(item_lst))):
            item1 = item_lst[item_index1]
            matrix_index1 = index_map[item1]
            for item_index2 in range(item_index1, len(item_lst)):
                item2 = item_lst[item_index2]
                matrix_index2 = index_map[item2]
                
                value = len([cluster for cluster in cluster_lst if item1 in cluster and item2 in cluster]) / partition_count
                matrix[matrix_index1, matrix_index2] = value
                matrix[matrix_index2, matrix_index1] = value
                
        if matrix.max() == 0:
            raise Exception("DEBUG THING")
                
        return CoAssosiationMatrix(matrix, index_map, partition_count)
    
    
    
    def GetEntry(self, item1: object, item2: object) -> float:
        Assert.assert_key_exists(item1, self.index_map)
        Assert.assert_key_exists(item2, self.index_map)
        i1, i2 = self.index_map[item1], self.index_map[item2]
        return self.matrix[i1, i2]
    
    def __getitem__(self, tuple: Tuple[object, object]) -> float:
        item1, item2 = tuple
        return self.GetEntry(item1, item2)