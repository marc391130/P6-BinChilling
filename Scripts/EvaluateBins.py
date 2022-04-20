import sys
from typing import List, Dict, Tuple
from tqdm import tqdm
from CheckMFilter import BinDto, CheckMFilter

sys.path.insert(1, '../ACE')
from Cluster import Cluster
from Domain import ContigData
from ContigReader import ContigReader
from BinEvaluator import BinEvaluator, ClusterReader


if __name__ == '__main__': 
    if len(sys.argv) != 8:
        print("arguments need to be:\n", \
            "1: Fasta filepath\n", \
            "2: SCG_Filepath\n", \
            "3: gene db file\n",\
            "4: Numpy_filepath\n", \
            "5: Cluster_filepath\n", \
            "6: depth_filepath\n",\
            "7: Output path")
    else:
        reader = ContigReader(fasta_file=sys.argv[1], depth_file=sys.argv[6], \
            SCG_filepath=sys.argv[2], SCG_db_path=sys.argv[3], numpy_file=sys.argv[4])
        cluster_reader = ClusterReader(file_path=sys.argv[5], contig_reader=reader)
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
        
        checkMfilter = CheckMFilter(None, sys.argv[7], lambda x: True)
        checkMfilter.write_output(result_dto_lst)
