from matplotlib.pyplot import axis
import numpy as np
from sklearn import cluster
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from typing import Tuple, Dict, List
import Assertions as Assert 

class MemberSimularityMatrix:
    def __init__(self, matrix: np.matrix, cluster_index_map: Dict[Cluster, int], item_index_map: Dict[Cluster, int]) -> None:
        self.matrix = matrix
        self.cluster_index_map = cluster_index_map
        self.item_index_map = item_index_map
    
    def shape(self) -> Tuple[int, int]:
        return (len(self.cluster_index_map), len(self.item_index_map))
    
    def __getitem__(self, tuple: Tuple[Cluster, object]):
        cluster, item = tuple
        return self.GetEntry(cluster, item)
    
    def Cluster_Mean(self, cluster: Cluster) -> float:
        Assert.assert_key_exists(cluster, self.cluster_index_map)
        c_index, count, value = self.cluster_index_map[cluster], len(self.cluster_index_map), 0
        for item, e_index in self.item_index_map.items():
            value += self.matrix[e_index, c_index]
        return value / count
    
    def Cluster_Max(self, cluster: Cluster) -> float:
        c_index, max_value = self.cluster_index_map[cluster], 0
        for item, e_index in self.item_index_map.items():
            max_value = max(self.matrix[e_index, c_index], max_value)
        return max_value
    
    def GetEntry(self, cluster: Cluster, item: object) -> float:
        Assert.assert_key_exists(cluster, self.cluster_index_map)
        Assert.assert_key_exists(item, self.item_index_map)
        c_index, e_index = self.cluster_index_map[cluster], self.item_index_map[item]
        return self.matrix[c_index, e_index]

class MemberMatrix:
    def __init__(self, clusters: List[Cluster], items: List) -> None:
        self.item_index_map = {items[i]: i for i in range(len(items))}
        self.initialize(clusters)
        
    
    #Overrides the internal matrix using only the clusters in the provided list
    def initialize(self, clusters: List[Cluster]):
        self.cluster_index_map = {clusters[i]: i for i in range(len(clusters))}
        self.membership_matrix : np.matrix = np.full(shape=self.shape(), fill_value=np.nan )
        self.__build_matrix__()
        
    def __build_matrix__(self) -> None:
        for cluster, c_index in self.cluster_index_map.items():
            membership_map = cluster.calc_all_membership()
            self.__calc_column__(c_index, membership_map)
    
    def __calc_column__(self, cluster_index: int,  membership_values: Dict[object, int]):
        for item, e_index in self.item_index_map.items():
                self.membership_matrix[cluster_index, e_index] = membership_values[item] if item in membership_values else 0
        
    def __len__(self) -> int:
        return len(self.cluster_index_map)
    
    def add_cluster(self, cluster: Cluster) -> None:
        item_map = cluster.calc_all_membership()
        column = [[item_map[item]] for item, i in self.item_index_map.items()]
        self.cluster_index_map[cluster] = len(self.cluster_index_map)
        
        self.membership_matrix = np.append(self.membership_matrix, [column], axis=1)
    
    def BuildMembershipSimilarityMatrix(self, clusters: set[Cluster]) -> MemberSimularityMatrix:
        cluster_lst = list(self.cluster_index_map.keys())
        merged_cluster_map = {cluster_lst[i]: i for i in range(len(cluster_lst)) if cluster_lst[i] in clusters}
        shape = (len(self.item_index_map), len(merged_cluster_map)) 
        matrix = np.full(shape=shape, fill_value=np.nan)
        
        def find_correct_max_membership(item_index) -> float:
            max_value = 0
            for cluster, c_index in self.cluster_index_map.items():
                if cluster in merged_cluster_map:
                    max_value = max(self.membership_matrix[item_index, c_index], max_value)
            return max_value
        
        for item, e_index in self.item_index_map.items():
            max_value = find_correct_max_membership(e_index)
            for cluster, c_index in merged_cluster_map.items():
                old_index = self.item_index_map[cluster]
                membervalue = self.membership_matrix[old_index, e_index]
                matrix[e_index, c_index] = membervalue / max_value
        return MemberSimularityMatrix(matrix, merged_cluster_map, self.item_index_map)        
        
    
    def shape(self) -> Tuple[int, int]:
        return (len(self.cluster_index_map), len(self.item_index_map))
    


        
                