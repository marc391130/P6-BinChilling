# Auther: Narucs, Dak, Spand, & Trudsliv
# Created: 23-02-2022
# Intended purpose: To get data from fasta files

from random import randrange, random
from AdaptiveEnsembler import AdaptiveClusterEnsembler
from Cluster import Cluster, Partition, PartitionSet
from tqdm import tqdm

ensembler = AdaptiveClusterEnsembler()

def dummy_data() -> PartitionSet:
    data = [random() for x in tqdm(range(1000))]

    partition_set = PartitionSet()
    for x in range(4):
        partition_set.append(Partition(data))

    for x in tqdm(data):
        rnd = randrange(0, 16)

        partition = int(rnd / 4)
        cluster = (int(rnd) % 4)

        partition_set[partition].add(str(cluster), x)
    return partition_set

ensembler.ensemble(dummy_data())



