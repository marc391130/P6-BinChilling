import numpy as np
from tqdm import tqdm
from Cluster import Cluster, PartitionSet
from typing import Tuple, Dict, List

class BitMatrix:
    def __init__(self, gamma: PartitionSet):
        self.matrix = None
        self.value_dct = None
        self.cluster_dct = None 
        self.bit_matrix_transform(gamma)

        pass

    def bit_matrix_transform(self, gamma: PartitionSet) -> None:
        clusters = gamma.get_all_clusters()
        self.value_dct = gamma.get_all_elements()
        
        self.cluster_dct = self.__setup_cluster_dct__(clusters)
        self.matrix = self.__setup_bit_matrix__(self.cluster_dct, self.value_dct)

        self.__update_matrix__()


    def add_cluster_to_bit_matrix(self, cluster: Cluster) -> None:
        self.cluster_dct[cluster] = len(self.cluster_dct)
        self.matrix = self.__setup_bit_matrix__(self.cluster_dct, self.value_dct)

        self.__update_matrix__()

    def __update_matrix__(self) -> None:
        for cluster, cluster_idx in self.cluster_dct.items():
            for item in cluster:
                self.matrix[self.value_dct[item], cluster_idx] = 1

    def calc_membership_in_cluster(self, item, cluster: Cluster) -> int:
        result = 0
        for leaf in cluster.__get_leaf_clusters__():
            result += self.get_entry(item, leaf)
        return result
        
    def membership_similarity_measure(self, item, cluster: Cluster, cluster_dct: Dict[Cluster, int]) -> float:
        # all_merged_cluster = dict(filter(lambda x: x[1] == MERGED_CLUSTER, cluster_dct.items())) if use_filter else cluster_dct
        # MIGHT HAVE FUCKED UP ^^^ BUT IT works maybe. Maybe filter is required, however it fixes bugs to not use.
        max_value = 0

        for cluster__, partition_idx in cluster_dct.items():
            max_value = max(self.calc_membership_in_cluster(item, cluster__), max_value)

        if max_value == 0:
            #return 0
            raise Exception(f"The item {str(item)} does not exist in the dataset")

        return self.calc_membership_in_cluster(item, cluster) / max_value

    def get_entry(self, item, cluster: Cluster) -> int:
        return self.matrix[self.value_dct[item], self.cluster_dct[cluster]]

    def __setup_cluster_dct__(self, clusters: List[Cluster]) -> Dict[Cluster, int]:
        result = {}
        for cluster_idx in range(len(clusters)):
            result[clusters[cluster_idx]] = cluster_idx
        return result


    def __setup_bit_matrix__(self, clusters, elements) -> np.matrix:
        column_nr = len(clusters) + 1
        row_nr = len(elements)

        matrix = np.empty_like(0, shape = (row_nr, column_nr))
        matrix.fill(0)
        return matrix

