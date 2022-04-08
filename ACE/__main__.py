# Auther: Narucs, Dak, Spand, & Trudsliv
# Created: 23-02-2022
# Intended purpose: To get data from fasta files

from random import randrange, random, seed

from AdaptiveEnsembler import AdaptiveClusterEnsembler, Ensembler, target_bin_3_4th_count_estimator
from Cluster import Cluster, Partition, PartitionSet, Contig
from tqdm import tqdm
from PartitionSetReader import PartitionSetReader
from ContigReader import ContigReader
import Constants
import time
import Constants as Const 


# import numpy
# arr = numpy.zeros((3,3))
# app = numpy.ones(shape=(3,1))
# print("> arr\n", arr)
# print("> app\n", app)
# print(numpy.append(arr , app, axis=1 ))
# raise Exception()

def print_partition(file_path: str, parititon: Partition):
    with open(file_path, 'w') as file:
        cluster_lst = list(parititon.values())
        for cluster_idx in range(len(cluster_lst)):
            for item in cluster_lst[cluster_idx]:
                file.write(f"{cluster_idx+1}\t{item.name}\n")

if __name__ == '__main__':

    # ensembler = AdaptiveClusterEnsembler( \
    #     initial_alpha1_thredshold=0.6, \
    #     initial_delta_aplha=0.1, \
    #     alpha1_min=0.9, \
    #     alpha2=0.6)
    ensembler = AdaptiveClusterEnsembler( \
        initial_alpha1_thredshold=0.8, \
        initial_delta_aplha=0.02, \
        alpha1_min=0.8, \
        alpha2=0.9,
        taget_clusters_est=target_bin_3_4th_count_estimator)

    seed(2) # most random number used for seed, chosen by committee

    use_real_data = True
    candidate_clusters = None
    if use_real_data:
        contigReader = ContigReader(Const.FASTA_FILEPATH, '../Dataset/oral_abundance.npz', Const.SCG_FILEPATH, "../Dataset/contigs_numpy.npy")
        contigReader.read_file(contigReader.fasta_file)
        partitionSetReader = PartitionSetReader("../Dataset/ClusterData/", contigReader, lambda x: x.endswith(".tsv"))
        partition_set = partitionSetReader.read_file()

        candidate_clusters = ensembler.ensemble(partition_set)
    else:
        def dummy_data(partitions = 2, clusters_per_partition = 100, elements_in_data=1000) -> PartitionSet:
            data = [Contig(str(random())) for x in tqdm(range(elements_in_data))]

            partition_set = PartitionSet(data)
            for x in range(partitions):
                partition_set.create_partition()

            for partition in range(len(partition_set)):
                cluster_max = max(int(randrange(int(clusters_per_partition * 0.1), clusters_per_partition)), 2)
                
                for x in data:
                    cluster = int(randrange(0,clusters_per_partition))
                    partition_set[partition].add(str(cluster), x)
            return partition_set

        candidate_clusters = ensembler.ensemble(dummy_data(partitions=10, clusters_per_partition=2, elements_in_data=10))

    print_partition("../Dataset/ACE_output.tsv", candidate_clusters)
    # ensembler.print_to_file("../Dataset/ACE_output.tsv", candidate_clusters)


