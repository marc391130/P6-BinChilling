from os import listdir
from os.path import join
from Cluster import Cluster, Partition, PartitionSet
from typing import Callable, Tuple
from ContigReader import ContigReader
import Constants as CONSTANT
from ContigData import ContigData
from tqdm import tqdm


class PartitionSetReader:
    def __init__(self, numpy_path: str,  folder_path: str, file_predicate: Callable[[str], bool]) -> None:
        self.folder_path = folder_path
        self.file_predicate = file_predicate
        self.contig_names = ContigReader().read_contig_names(CONSTANT.FILEPATH)
        self.contigs = ContigReader().read_file_fast(CONSTANT.FILEPATH, numpy_path)

    def read_partisionSet(self) -> PartitionSet:
        return self.read_file()

    def read_file(self, show_warnings: bool = True) -> PartitionSet:
        files = [f for f in listdir(self.folder_path) if self.file_predicate(f)]
        
        partition_set = PartitionSet[ContigData](list(self.contigs.values()))
        
        warning_clusters = []
        print("reading partitions files...")
        for file in tqdm(files):
            partition = partition_set.create_partition() #dont have to add it to partition after this
            self.__parse_partition__(file, partition)
            
            if len(partition) <= 1:
                warning_clusters.append(f"partition in file { file } only contains 1 cluster")
        
        if show_warnings:
            [print("> " +w) for w in warning_clusters]
        return partition_set

    def __parse_partition__(self, filename: str, partition: Partition) -> None:
        with open(join(self.folder_path, filename), 'r') as f:
            for line in f.readlines():
                cluster_name, edge_name = self.__parse_cluster_line__(line) 
                partition.add(cluster_name, self.contigs[edge_name])

                
    #returns tuple of cluster_name, edge_name
    def __parse_cluster_line__(self, line: str) -> Tuple[str, str]:
        split_line = line.split('\t')
        return (split_line[0], split_line[1].replace('\n', ''))
