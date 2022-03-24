# Auther: Narucs, Dak, Spand, & Trudsliv
# Created: 23-02-2022
# Intended purpose: To get data from fasta files

from random import randrange, random, seed
from AdaptiveEnsembler import AdaptiveClusterEnsembler
from Cluster import Cluster, Partition, PartitionSet, Contig
from tqdm import tqdm
from PartitionSetReader import PartitionSetReader
from ContigReader import ContigReader
import Constants
import time

if __name__ == '__main__':

    ensembler = AdaptiveClusterEnsembler( \
        initial_alpha1_thredshold=0.6, \
        initial_delta_aplha=0.1, \
        alpha1_min=0.20, \
        alpha2=0.1)

    seed(2)

    use_real_data = True
    candidate_clusters = None
    if use_real_data:
        partitionSetReader = PartitionSetReader("../Dataset/contigs_numpy.npy", "../Dataset/ClusterData/", lambda x: x.endswith(".tsv"))
        partition_set = partitionSetReader.read_file()

        candidate_clusters = ensembler.ensemble(partition_set)
    else:
        def dummy_data(partitions = 2, clusters_per_partition = 100, elements_in_data=1000) -> PartitionSet:
            data = [Contig(str(random())) for x in tqdm(range(elements_in_data))]

            partition_set = PartitionSet(data)
            for x in range(partitions):
                partition_set.create_partition()

            for partition in range(len(partition_set)):
                for x in data:
                    cluster = int(randrange(0,clusters_per_partition))
                    partition_set[partition].add(str(cluster), x)
            return partition_set

        candidate_clusters = ensembler.ensemble(dummy_data(partitions=2, clusters_per_partition=100, elements_in_data=1000))

    ensembler.print_to_file("../Dataset/ACE_output.tsv", candidate_clusters)


