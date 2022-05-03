import itertools
import sys
from typing import Iterator, List, Dict, Tuple, TypeVar, Generic
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
            SCGs_count = self.__calculate_number_of_SCGs__(cluster, skip_item, None)
            result[cluster] = (self.__calculate_completeness__(SCGs_count), self.__calculate_contamination__(SCGs_count))
            #print((self.__calculate_completeness__(cluster)) - (0.5 * self.__calculate_contamination__(cluster)) - (0.5 * self.__calculate_megabin_penalty__(cluster)))

        return result

    def evaluate_cluster(self, cluster: Cluster, skip_item: ContigData = None) -> Tuple[float, float, float]:
        SCGs_count = self.__calculate_number_of_SCGs__(cluster, skip_item, None)
        return (self.__calculate_completeness__(SCGs_count), self.__calculate_contamination__(SCGs_count), self.__calculate_megabin_penalty__(SCGs_count))
        
    def evaluate_item_lst(self, item_set: set):
        SCGs_count = self.__calculate_number_of_SCGs__(item_set, None, None)
        return (self.__calculate_completeness__(SCGs_count), self.__calculate_contamination__(SCGs_count), self.__calculate_megabin_penalty__(SCGs_count)) 
        
    

    def calculate_score(self, cluster: Cluster, skip_item: ContigData = None, include_item: ContigData = None) -> float:
        scg_count = self.__calculate_number_of_SCGs__(cluster, skip_item, include_item)
        return self.__score_from_SCG_count__(scg_count)
    
    
    def __score_from_SCG_count__(self, SCG_count: Dict[str, int]) -> float:
        completeness, contamination, megabin =\
                (self.__calculate_completeness__(SCG_count),\
                self.__calculate_contamination__(SCG_count),\
                self.__calculate_megabin_penalty__ (SCG_count))
                
        return self.score_calc(completeness, contamination, megabin)

    def score_calc(self, completeness: float, contamination: float, megabin_pen: float):
        return completeness - max(contamination**2, contamination) - (0.4*megabin_pen)

    def calculate_item_score(self, cluster: Cluster, extra_item: ContigData or None = None) -> Dict[ContigData, float]:
        if len(cluster) == 0: return {}
        if len(cluster) == 1: 
            return { item: self.calculate_score(cluster, include_item=extra_item) for item in self.__chain_cluster__(cluster, extra_item) }
        SCGs = self.__calculate_number_of_SCGs__(cluster, None, extra_item)
        total_value = self.__score_from_SCG_count__(SCGs)

        result = {item: total_value-  self.__score_from_SCG_count__(self.__remove_item_from_SCG_Count__(item, SCGs))\
            for item in self.__chain_cluster__(cluster, extra_item)}
        return result
    

    
    def __remove_item_from_SCG_Count__(self, item: ContigData, scgs: Dict[str, int]) -> Dict[int, str]:
        return { scg: ((count - 1) if scg in item.SCG_genes else count )\
            for scg, count in scgs.items() if not (count <= 1 and scg in item.SCG_genes)  }

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

    def __calculate_completeness__(self, SCG_count: Dict[str, int]) -> float:
        
        counter = len(SCG_count)
        divisor = len(self.all_SCGs)

        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_contamination__(self, SCG_count: Dict[str, int]) -> float:
        # uniques = self.__calculate_unqiues__(cluster, skip_item, include_item)

        dSCG = [scg for scg, count in SCG_count.items() if count > 1]

        counter = len(dSCG)
        divisor = len(SCG_count) #This is the same as uniques, but does not require an additional n^2 call

        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_purity__(self, SCG_count: Dict[str, int]) -> float:
        #uniques = self.__calculate_unqiues__(cluster, skip_item, include_item)
        # SCGs = self.__calculate_number_of_SCGs__(cluster, skip_item, include_item)
        dSCG = [scg for scg, count in SCG_count.items() if count > 1]
        
        counter = len(SCG_count) #Same as uniques, without having to compute
        divisor = len(SCG_count) + len(dSCG)
        
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_megabin_penalty__(self, SCG_count: Dict[str, int]) -> float:
        # SCGs = self.__calculate_number_of_SCGs__(cluster, skip_item, include_item)
        #uniques = self.__calculate_unqiues__(cluster, skip_item, include_item)
        nr_SCGs = sum(SCG_count.values()) - len(SCG_count)

        return (nr_SCGs / len(self.all_SCGs))*100 if len(self.all_SCGs) != 0 else 0

    def __calculate_unqiues__(self, cluster: Cluster[ContigData], skip_item: ContigData = None, include_item: ContigData = None) -> set:
        uniques = set()

        for item in self.__chain_cluster__(cluster, include_item):
            if skip_item is not None and item is skip_item: continue
            uniques.update(item.SCG_genes)

        return uniques
    
    def __calculate_number_of_SCGs__(self, cluster: Cluster[ContigData], skip_item: ContigData = None, include_item: ContigData = None) -> Dict[str, int]:
        SCGs = {}
        
        for item in self.__chain_cluster__(cluster, include_item):
            if skip_item is not None and item is skip_item: continue
            for scg in item.SCG_genes:
                SCGs[scg] = SCGs[scg] + 1 if scg in SCGs else 1

        return SCGs
    
    def __chain_cluster__(self, cluster: Cluster, include_item: ContigData or None = None) -> Iterator[ContigData]:
        return cluster.__iter__() if include_item is None or include_item in cluster\
            else itertools.chain(cluster.__iter__(), [include_item])

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