from typing import Dict, List, Tuple
from Cluster import Cluster, Contig, Partition, PartitionSet
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

    def generate_consensus_clusters(self, gamma: PartitionSet, bit_matrix: np.matrix) -> List[Cluster]:

        dct_info = {}
        merged = []

        for partition_idx in range(len(gamma)):
            for key, value in gamma[partition_idx].items():
                dct_info[value] = partition_idx

        done_dct = {}

        for key1, value1 in dct_info.items():
            if key1 in done_dct:
                continue
            for key2, value2 in dct_info.items():
                if key1 is key2 or value1 == value2 or key2 in done_dct:
                    continue
            
                similarity = gamma.__similarity_measure_cluster__(key1, key2)

                if similarity >= self.alpha1_thredshold:
                    merged.append(key1.merge(key2))
                    done_dct[key1] = True
                    done_dct[key2] = True

        similarities_found = 1
        while similarities_found != 0:
            similarities_found = 0
            done_dct = {}
            for key, value in dct_info.items():
                for cluster in merged:
                    similarity = gamma.__similarity_measure_cluster__(key, cluster)
                    if similarity >= self.alpha1_thredshold:
                        done_dct[key] = 0
                        merged.append(key.merge(cluster))
                        merged.remove(cluster)
                        similarities_found += 1

            for key, value in done_dct.items():
                del dct_info[key]
            print(similarities_found, len(dct_info), len(merged))
        
        print(f"Resulted in lambda merged clusters: {len(merged)}")

        return merged                


    def __membership_similarity__(self, bit_matrix: np.matrix, item, cluster: Cluster, merged_lst: List[Cluster]) -> None:
        pass


    def ensemble(self, gamma: PartitionSet[Contig]) -> None:
        data = self.bit_matrix_transform(gamma)
        print(data, data.shape)
        merged = self.generate_consensus_clusters(gamma, data)
        print(merged)

