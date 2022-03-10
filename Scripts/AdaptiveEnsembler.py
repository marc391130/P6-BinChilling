from typing import Dict, List, Tuple
from Cluster import Contig, Partition, PartitionSet
import numpy as np
from tqdm import tqdm

class Ensembler:
    def ensemble(self, data) -> None:
        pass

class AdaptiveClusterEnsembler(Ensembler):
    def __init__(self, initial_alpha1_thredshold = 0.8):
        self.alpha1_thredshold = initial_alpha1_thredshold
        pass

    def bit_matrix_transform(self, gamma: PartitionSet) -> np.matrix:
        matrix = self.__setup_bit_matrix(gamma)
        
        values_idx = gamma.get_all_elements()
        clusters = gamma.get_all_clusters()
        for cluster_idx in tqdm(range(len(clusters))):
            for item in clusters[cluster_idx]:
                matrix[values_idx[item], cluster_idx] = True
        
        return matrix

    def __setup_bit_matrix(self, gamma: PartitionSet) -> np.matrix:
        column_nr = len(gamma.get_all_clusters())
        row_nr = len(gamma.get_all_elements())

        matrix = np.empty_like(False, shape = (row_nr, column_nr))
        matrix.fill(False)
        return matrix

    def ensemble(self, gamma: PartitionSet[Contig]) -> None:
        data = self.bit_matrix_transform(gamma)
        print(data, data.shape)

