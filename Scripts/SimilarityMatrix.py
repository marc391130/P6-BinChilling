import numpy as np
from torch import float32
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from typing import Tuple, Dict, List

class SimilarityMatrix(Dict[Cluster, Dict[Cluster]]):
    def __init__(self) -> None:
        self.cluster_count = 0
        self.matrix = None
    
    
    def initialize(self, gamma: PartitionSet) -> None:
        self.cluster_count = len(gamma.get_all_clusters())
        self.matrix = np.matrix(shape=(self.cluster_count,self.cluster_count), dtype=float32)