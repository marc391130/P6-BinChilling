# Auther: Narucs, Dak, Spand, & Trudsliv
# Created: 23-02-2022
# Intended purpose: To get data from fasta files

from random import randrange, random, seed
from AdaptiveEnsembler import AdaptiveClusterEnsembler
from Cluster import Cluster, Partition, PartitionSet
from tqdm import tqdm
from PartitionSetReader import PartitionSetReader
from ContigReader import ContigReader
import Constants

ensembler = AdaptiveClusterEnsembler(0.6, 0.1, 0.4, 0.6)

seed(1)

partitionSetReader = PartitionSetReader("../Dataset/ClusterData/", lambda x: x.endswith(".tsv"))
partition_set = partitionSetReader.read_file()


ensembler.ensemble(partition_set)
# def dummy_data(partitions = 4, clusters_per_partition = 10, elements_in_data=1000) -> PartitionSet:
#     data = [random() for x in tqdm(range(elements_in_data))]

#     partition_set = PartitionSet()
#     for x in range(partitions):
#         partition_set.append(Partition(data))

#     for partition in range(len(partition_set)):
#         for x in data:
#             cluster = int(randrange(0,clusters_per_partition))
#             partition_set[partition].add(str(cluster), x)
#     return partition_set

#ensembler.ensemble(dummy_data(partitions=10, clusters_per_partition=10, elements_in_data=100))



