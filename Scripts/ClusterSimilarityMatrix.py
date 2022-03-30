import numpy as np
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from typing import Tuple, Dict, List
import Assertions as Assert 

class ClusterSimilarityMatrix():
    def __init__(self, matrix: np.matrix, index_map: Dict[Cluster, int]) -> None:
        print(">SIZE", np.shape(matrix))
        self.matrix = matrix
        self.index_map = index_map
    
    def __len__(self) -> int:
        return len(self.index_map)

        
    def __getitem__(self, tuple: Tuple[Cluster, Cluster]) -> int or np.nan:
        c1, c2 = tuple
        return self.getEntry(c1, c2)
    
    def __setitem__(self, object: Cluster) -> None:
        raise Exception("SimilarityMatrix is not allowed to be updated")
    
    
    def __expand__(self, cluster: Cluster, values: Dict[int, float or np.nan]) -> None:
        Assert.assert_key_exists(cluster, self.index_map)
        new_index = self.index_map[cluster]
        size = len(values)
        
        new_row = np.empty_like(np.nan, shape=(1, size-1))
        new_column = np.empty_like(np.nan, shape=(size, 1))
        
        self.matrix = np.append(self.matrix, new_row,    axis=0 )       
        self.matrix = np.append(self.matrix, new_column, axis=1 )
        
        for index, value in values.items():
            self.matrix[new_index, index] = value
            self.matrix[index, new_index] = value
        return
    
    def ToDict(self, for_cluster: Cluster) -> Dict[Cluster, int or np.nan]:
        Assert.assert_key_exists(for_cluster, self.index_map)
        index = self.index_map[for_cluster]
                
        return { cluster: self.getByIndex(index, other_index) for cluster, other_index in self.index_map.items()  }     
    
    
    def getEntry(self, cluster1 : Cluster, cluster2: Cluster) -> int or np.nan:
        id1, id2 = self.index_map[cluster1], self.index_map[cluster2]
        return self.getByIndex(id1, id2)
    
    def getByIndex(self, index1: int, index2: int) -> int or np.nan:
        if index1 == index2:
            return np.nan
        return self.matrix[index1, index2]
    
    def max_similarity(self, clusters: List[Cluster] = None) -> float:
        if clusters is None:
            return self.matrix.max()
        
        max_value = np.nan
        for cl_index1 in range(len(clusters)):
            index1 = self.index_map[clusters[cl_index1]]
            for cl_index2 in range(cl_index1, len(clusters)):
                index2 = self.index_map[clusters[cl_index2]]
                value = self.getByIndex(index1, index2)
                if np.isnan(value):
                    continue
                elif np.isnan(max_value):
                    max_value = value 
                else:
                    max_value = max(value, max_value)
        return max_value
            
        