import sys
from typing import List, Dict, Tuple
from tqdm import tqdm

from Cluster import Cluster
from Domain import ContigData
from ContigReader import ContigReader


class BinEvaluator:
    def __init__(self, all_SCGs: set) -> None:
        self.all_SCGs = all_SCGs

    def score(self, cluster_lst: List[Cluster]) -> Dict[Cluster, float]:
        return { cluster: self.calculate_score(cluster) for cluster in cluster_lst }


    def evaluate(self, cluster_lst: List[Cluster], skip_item: ContigData = None) -> Dict[Cluster, Tuple[float, float]]:
        result = {}

        for cluster in cluster_lst:
            result[cluster] = (self.__calculate_completeness__(cluster, skip_item), self.__calculate_contamination__(cluster, skip_item))
            #print((self.__calculate_completeness__(cluster)) - (0.5 * self.__calculate_contamination__(cluster)) - (0.5 * self.__calculate_megabin_penalty__(cluster)))

        return result

    def calculate_score(self, cluster: Cluster, skip_item: ContigData = None) -> float:
        completeness, contamination, megabin =\
                (self.__calculate_completeness__(cluster, skip_item),\
                self.__calculate_contamination__(cluster, skip_item),\
                #0,\
                self.__calculate_megabin_penalty__ (cluster, skip_item))
                
        return completeness - (0.5*contamination) - (megabin*0.5)

    def __calculate_sight__(self, completeness, contamination) -> str:
        if completeness > 90 and contamination < 5:
            return 'near'
        if contamination > 5 and contamination <= 10 and\
            completeness < 90 and completeness >= 70:
            return 'substantial'
        if contamination > 10 and contamination <= 15 and\
            completeness < 70 and completeness >= 50:
            return 'moderate'
        if contamination > 15 and completeness < 50:
            return 'partial'
        if (completeness - contamination) > 0: return 'bad'
        return 'zero'

    def __calculate_completeness__(self, cluster: Cluster[ContigData], skip_item: ContigData = None) -> float:
        uniques = self.__calculate_unqiues__(cluster, skip_item)
        
        counter = len(uniques)
        divisor = len(self.all_SCGs)

        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_contamination__(self, cluster: Cluster[ContigData], skip_item: ContigData = None) -> float:
        SCGs = self.__calculate_number_of_SCGs__(cluster, skip_item)
        uniques = self.__calculate_unqiues__(cluster, skip_item)

        dSCG = [scg for scg, count in SCGs.items() if count > 1]

        counter = len(dSCG)
        divisor = len(uniques)

        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_purity__(self, cluster: Cluster[ContigData], skip_item: ContigData = None) -> float:
        uniques = self.__calculate_unqiues__(cluster, skip_item)
        SCGs = self.__calculate_number_of_SCGs__(cluster, skip_item)
        dSCG = [scg for scg, count in SCGs.items() if count > 1]
        
        counter = len(uniques)
        divisor = len(uniques) + len(dSCG)
        
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_megabin_penalty__(self, cluster: Cluster[ContigData], skip_item: ContigData = None) -> float:
        SCGs = self.__calculate_number_of_SCGs__(cluster, skip_item)
        uniques = self.__calculate_unqiues__(cluster, skip_item)
        nr_SCGs = sum(list(SCGs.values())) - len(uniques)

        return (nr_SCGs / len(self.all_SCGs))*100 if len(self.all_SCGs) != 0 else 0

    def __calculate_unqiues__(self, cluster: Cluster[ContigData], skip_item: ContigData = None) -> set:
        uniques = set()

        for item in cluster:
            if skip_item is not None and item is skip_item: continue
            for scg in item.SCG_genes:
                uniques.add(scg)

        return uniques
    
    def __calculate_number_of_SCGs__(self, cluster: Cluster[ContigData], skip_item: ContigData = None) -> Dict[str, int]:
        SCGs = {}

        for item in cluster:
            if skip_item is not None and item is skip_item: continue
            for scg in item.SCG_genes:
                SCGs[scg] = SCGs[scg] + 1 if scg in SCGs else 1

        return SCGs

class ClusterReader:
    def __init__(self, file_path: str, contig_reader: ContigReader) -> None:
        self.contigReader = contig_reader
        self.clusters = self.__read_clusters__(file_path)

    def __read_clusters__(self, file_path: str) -> List[Cluster]:
        result_clusters = []
        cluster_data_map = {}
        with open(file_path, 'r') as f:
            for line in tqdm(f.readlines()):
                split_line = line.split('\t')
                cluster_idx, edge_name = split_line[0], split_line[1].replace('\n', '')
                cluster_data_map[cluster_idx] = cluster_data_map[cluster_idx] + [edge_name] if cluster_idx in cluster_data_map else [edge_name]
        
        contig_scg_dct = self.contigReader.read_file_fast(None, True)
        for cluster_idx, edge_lst in cluster_data_map.items():
            cluster = Cluster()
            for edge in edge_lst:
                r = contig_scg_dct.get(edge, None)
                if r is not None: cluster.add(r)
            result_clusters.append(cluster)
        return result_clusters