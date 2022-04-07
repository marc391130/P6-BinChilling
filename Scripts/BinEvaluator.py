import sys
from typing import List, Dict, Tuple
from tqdm import tqdm
from CheckMFilter import BinDto, CheckMFilter

sys.path.insert(1, "../ACE")

from Cluster import Cluster
from ContigData import ContigData
from ContigReader import ContigReader


class BinEvaluator:
    def __init__(self, all_SCGs: set) -> None:
        self.all_SCGs = all_SCGs

    def evaluate(self, cluster_lst: List[Cluster]) -> Dict[Cluster, Tuple[float, float]]:
        result = {}

        for cluster in cluster_lst:
            result[cluster] = (self.__calculate_completeness__(cluster), self.__calculate_contamination__(cluster))
            #print((self.__calculate_completeness__(cluster)) - (0.5 * self.__calculate_contamination__(cluster)) - (0.5 * self.__calculate_megabin_penalty__(cluster)))

        return result

    def __calculate_completeness__(self, cluster: Cluster[ContigData]) -> float:
        uniques = self.__calculate_unqiues__(cluster)
        
        counter = len(uniques)
        divisor = len(self.all_SCGs)

        return (counter / divisor) * 100 if divisor != 0 else 0

    def __calculate_contamination__(self, cluster: Cluster[ContigData]) -> float:
        SCGs = self.__calculate_number_of_SCGs__(cluster)
        uniques = self.__calculate_unqiues__(cluster)

        dSCG = [scg for scg, count in SCGs.items() if count > 1]

        counter = len(dSCG)
        divisor = len(uniques)

        return (counter / divisor) * 100 if divisor != 0 else 0

    def __calculate_megabin_penalty__(self, cluster: Cluster[ContigData]) -> float:
        SCGs = self.__calculate_number_of_SCGs__(cluster)
        uniques = self.__calculate_unqiues__(cluster)
        nr_SCGs = sum(list(SCGs.values())) - len(uniques)

        return nr_SCGs / len(self.all_SCGs) if len(self.all_SCGs) != 0 else 0

    def __calculate_unqiues__(self, cluster: Cluster) -> set:
        uniques = set()

        for item in cluster:
            for scg in item.SCG_genes:
                uniques.add(scg)

        return uniques
    
    def __calculate_number_of_SCGs__(self, cluster: Cluster[ContigData]) -> Dict[str, int]:
        SCGs = {}

        for item in cluster:
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
                cluster.add(contig_scg_dct[edge])
            result_clusters.append(cluster)
        return result_clusters

                
if __name__ == '__main__': 
    if len(sys.argv) != 5:
        print("arguments need to be:\n", \
            "1: SCG_Filepath\n", \
            "2: Numpy_filepath\n", \
            "3: Cluster_filepath", \
            "4: Output path")
    else:
        reader = ContigReader("", "", sys.argv[1], sys.argv[2])
        cluster_reader = ClusterReader(sys.argv[3], reader)
        clusters = cluster_reader.clusters

        all_scgs = reader.read_total_SCGs_set() # 
        evaluator = BinEvaluator(all_scgs)
        data = evaluator.evaluate(clusters)

        result_dto_lst = []
        data_lst = list(data.items())
        for data_idx in range(len(data_lst)):
            cluster_idx, tuple = data_lst[data_idx]
            completeness, contamination = tuple

            dto = BinDto(f"bin_{data_idx}", contamination, completeness)
            result_dto_lst.append(dto)
        
        checkMfilter = CheckMFilter(None, sys.argv[4], lambda x: True)
        checkMfilter.write_output(result_dto_lst)


    
    