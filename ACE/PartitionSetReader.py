from os import listdir
from os.path import join
from Cluster import Cluster, Partition, PartitionSet
from typing import Callable, Tuple, Dict, List
from ContigReader import ContigReader
import Constants as CONSTANT
from ContigData import ContigData
from tqdm import tqdm


class PartitionSetReader:
    def __init__(self, cluster_folder_path: str, contig_reader: ContigReader, file_predicate: Callable[[str], bool] = lambda x: True) -> None:
        self.folder_path = cluster_folder_path
        self.file_predicate = file_predicate
        self.contig_reader = contig_reader

    def read_partisionSet(self) -> PartitionSet:
        return self.read_file()

    def read_file(self, show_warnings: bool = True) -> PartitionSet:
        contig_dct = self.contig_reader.read_file_fast(None, True)
        files = [f for f in listdir(self.folder_path) if self.file_predicate(f)]
        
        partition_set = PartitionSet[ContigData](list(contig_dct.values()))
        
        warning_clusters = []
        print("reading partitions files...")
        for file in tqdm(files):
            partition = partition_set.create_partition() #dont have to add it to partition after this
            self.__parse_partition__(file, partition, contig_dct)
            
            if len(partition) <= 1:
                warning_clusters.append(f"partition in file { file } only contains 1 cluster")
        
        if show_warnings:
            [print("> " +w) for w in warning_clusters]
        return partition_set

    def __parse_partition__(self, filename: str, partition: Partition, contigs: Dict[str, ContigData]) -> None:
        with open(join(self.folder_path, filename), 'r') as f:
            for line in f.readlines():
                cluster_name, edge_name = self.__parse_cluster_line__(line) 
                partition.add(cluster_name, contigs[edge_name])

                
    #returns tuple of cluster_name, edge_name
    def __parse_cluster_line__(self, line: str) -> Tuple[str, str]:
        split_line = line.split('\t')
        return (split_line[0], split_line[1].replace('\n', ''))
