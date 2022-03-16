from os import listdir
from os.path import join
from Cluster import Cluster, Partition, PartitionSet
from typing import Callable
from ContigReader import ContigReader
import Constants as CONSTANT

class PartitionSetReader:
    def __init__(self, folder_path: str, file_predicate: Callable[[str], bool]) -> None:
        self.folder_path = folder_path
        self.file_predicate = file_predicate
        self.data = ContigReader().read_contig_names(CONSTANT.FILEPATH)

    def read_file(self) -> PartitionSet:
        files = [f for f in listdir(self.folder_path) if self.file_predicate(f)]
        
        partition_set = PartitionSet()
        for file in files:
            partition = Partition(self.data)
            with open(join(self.folder_path, file), 'r') as f:
                lines = f.readlines()
                for line in lines:
                    split_line = line.split('\t')
                    cluster, edge = split_line[0], split_line[1].replace('\n', '')
                    partition.add(cluster, edge)

            if len(partition) == 1:
                print(str(file))
            partition_set.append(partition)
        
        return partition_set

                
                
