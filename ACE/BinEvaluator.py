import itertools
import sys
from typing import Iterable, Iterator, List, Dict, Tuple, TypeVar, Generic
from tqdm import tqdm

from Cluster import Cluster
from Domain import ContigData
from ContigReader import ContigReader



class BinEvaluator:
    def __init__(self, all_SCGs: set, genome_size_range: Tuple[int, int]  ) -> None:
        self.all_SCGs = all_SCGs
        if genome_size_range[0] > genome_size_range[1]: raise Exception('Minimum genome size bigger than max.')
        self.min_genome_len, self.max_genome_len = genome_size_range
        self.avg_genome_size = (self.min_genome_len + self.max_genome_len) / 2
        
    def score(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> float:
        comp, conn, mp = self.evaluate(cluster, skip_item, include_item)
        return self.calc_score(comp, conn, mp)
    
    
    def score_size(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> float:
        total_size = self.__sum_size__(cluster, skip_item, include_item)
        if self.min_genome_len <= total_size <= self.max_genome_len: 
            return 1.0
        elif total_size < self.min_genome_len:
            return 1 - ( (1 / self.min_genome_len) * (total_size - self.min_genome_len)**2 )
        elif total_size > self.max_genome_len:
            return 1 - ( (1 / self.min_genome_len) * (self.min_genome_len - total_size)**2 )
        return 0
        
                
    
    def score_SCG(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> float:
        comp, conn, mp = self.evaluate(cluster, skip_item, include_item)
        return self.calc_score(comp, conn, mp)
    
    def evaluate_lst(self, clusters: List[Cluster]) -> Dict[Cluster, Tuple[float, float, float]]:
        return {cluster: self.evaluate(cluster) for cluster in clusters } 
    
    def evaluate(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> Tuple[float, float, float]:
        scg_count = self.__calculate_number_of_SCGs__(cluster, skip_item, include_item)
        return self.__evaluate_scg_count__(scg_count)
    
    
    def score_lst(self, cluster_lst: List[Cluster]) -> Dict[Cluster, float]:
        return {cluster: self.score(cluster) for cluster in cluster_lst}
    
    def score_item_lst(self, item_lst: Iterable[ContigData]) -> float:
        scg_count = self.__calculate_number_of_SCGs__(item_lst)
        return self.__score_scg_count__(scg_count)
    
    def calc_score(self, completeness: float, contamination: float, megabin_pen: float) -> float:
        return completeness - (contamination**2) - (0.4*megabin_pen)
        # return completeness - 0.5*contamination - 0.5*megabin_pen
    
    def score_items(self, cluster: Cluster, extra_item: ContigData = None) -> Dict[ContigData, float]:
        if len(cluster) == 0 and extra_item is None: return {}
        if len(cluster) == 1 and (extra_item is None or extra_item in cluster) or\
            (len(cluster) == 0 and extra_item is not None): 
            return { item: self.score(cluster, include_item=extra_item) for item in self.__chain_cluster__(cluster, extra_item) }
        SCGs = self.__calculate_number_of_SCGs__(cluster, None, extra_item)
        total_value = self.__score_scg_count__(SCGs)

        result = {item: total_value -  self.__score_scg_count__(self.__remove_item_from_SCG_count__(item, SCGs))\
            for item in self.__chain_cluster__(cluster, extra_item)}
        return result
    
    def calculate_item_score(self, cluster: Cluster, extra_item: ContigData or None = None) -> Dict[ContigData, float]:
        if len(cluster) == 0 and extra_item is None: return {}
        if len(cluster) == 1 and (extra_item is None or extra_item in cluster) or\
            (len(cluster) == 0 and extra_item is not None): 
            return { item: self.score(cluster, include_item=extra_item) for item in self.__chain_cluster__(cluster, extra_item) }
        SCGs = self.__calculate_number_of_SCGs__(cluster, None, extra_item)
        total_value = self.__score_scg_count__(SCGs)

        result = {item: total_value-  self.__score_scg_count__(self.__remove_item_from_SCG_count__(item, SCGs))\
            for item in self.__chain_cluster__(cluster, extra_item)}
        return result
    
    
    
    def __remove_item_from_SCG_count__(self, item: ContigData, scgs: Dict[str, int]) -> Dict[int, str]:
        return { scg: ((count - 1) if scg in item.SCG_genes else count )\
            for scg, count in scgs.items() if not (count <= 1 and scg in item.SCG_genes)  }
    
    #?______________ PRIVATE METHODS ______________
    def __score_scg_count__(self, scg_count: Dict[str, int]) -> float:
        comp, conn, mp = self.__evaluate_scg_count__(scg_count)
        return self.calc_score(comp, conn, mp)
    
    def __evaluate_scg_count__(self, scg_count: Dict[str, int]) -> Tuple[float,float,float]:
        return \
            self.__calculate_completeness__(scg_count), \
            self.__calculate_contamination__(scg_count),\
            self.__calculate_megabin_penalty__(scg_count)
    
    def __calculate_completeness__(self, SCG_count: Dict[str, int]) -> float:
        
        counter = len(SCG_count)
        divisor = len(self.all_SCGs)
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_contamination__(self, SCG_count: Dict[str, int]) -> float:
        dSCG = [scg for scg, count in SCG_count.items() if count > 1]

        counter = len(dSCG)
        divisor = len(SCG_count) #This is the same as uniques, but does not require an additional n^2 call
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_purity__(self, SCG_count: Dict[str, int]) -> float:
        dSCG = [scg for scg, count in SCG_count.items() if count > 1]
        
        counter = len(SCG_count) #Same as uniques, without having to compute
        divisor = len(SCG_count) + len(dSCG)
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_megabin_penalty__(self, SCG_count: Dict[str, int]) -> float:
        nr_SCGs = sum(SCG_count.values()) - len(SCG_count)
        return (nr_SCGs / len(self.all_SCGs))*100 if len(self.all_SCGs) != 0 else 0
    
    def __calculate_number_of_SCGs__(self, cluster: Iterable[ContigData], skip_item: ContigData = None, include_item: ContigData = None) -> Dict[str, int]:
        SCGs = {}
        for item in self.__chain_cluster__(cluster, include_item):
            if skip_item is not None and item is skip_item: continue
            for scg in item.SCG_genes: 
                SCGs[scg] = SCGs.get(scg, 0) + 1
        return SCGs
    
    def __sum_size__(self, cluster: Cluster, skip_item: ContigData or None = None, include_item: ContigData or None = None):
        return sum((x.contig_length for x in self.__chain_cluster__(cluster, include_item) if x is not skip_item))
    
    def __chain_cluster__(self, cluster: Iterable[ContigData], include_item: ContigData or None = None) -> Iterator[ContigData]:
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