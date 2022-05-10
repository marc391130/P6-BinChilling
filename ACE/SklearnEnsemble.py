import os
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from scipy.sparse import lil_matrix, csr_matrix, csc_matrix
import argparse
import numpy as np
from sklearn.model_selection import cross_val_score
from tqdm import tqdm

from Domain import ContigData
from Cluster import Cluster, Partition
from ContigReader import ContigReader
from PartitionSetReader import PartitionSetReader

def parser():
    parser = argparse.ArgumentParser(
        prog='SKLearn Ensemble test',
        description="""SKLearn Ensemble test""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s WRITE THIS LATER PLZ",
        add_help=True
        )

    io_args = parser.add_argument_group(title='Contig input (required)', description=None)
    io_args.add_argument('--fasta' ,'-f', metavar='', required=True,\
        dest='fasta', help='path to fasta file of contigs')
    io_args.add_argument('--cache' ,'-c', metavar='', required=False,\
        dest='cache', help='path to cache file (Give this feature as .npy format)')
    io_args.add_argument('--partitions' ,'-p', metavar='', required=True,\
        dest='partionFolder', help='path to folder of partitions')
    io_args.add_argument('--depth' ,'-d', metavar='', required=True,\
        dest='depth', help='path to depthfile')
    io_args.add_argument('--threads' ,'-t', metavar='', required=False,\
        dest='threads', help='path to depthfile', default=1, type=int)
    
    args = parser.parse_args()
    __assert_args__(args)
    return args

def __assert_args__(args):
    if not os.path.isdir(args.partionFolder):
        raise Exception(f"{args.partionFolder} is not a folder path!")

    partition_paths = [os.path.join(args.partionFolder, f) for f in os.listdir(args.partionFolder) if f.endswith('.tsv')]

    if not os.path.isfile(args.fasta):
        raise Exception(f"{args.fasta} is not a file!")

    for partition_path in partition_paths:
        if not os.path.isfile(partition_path):
            raise Exception(f"{partition_path} is not a file!")

    return None

def run(args):

    # Read Fasta!
    contig_reader = ContigReader(args.fasta, depth_file=args.depth, numpy_file=args.cache, enable_analyse_contig_comp=True, max_threads=args.threads)
    contigs = contig_reader.read_file_fast(args.cache)

    # Read Partitions!
    partition_reader = PartitionSetReader(args.partionFolder, contig_reader)
    partition_set = partition_reader.read_partisionSet()

    clusters = partition_set.get_all_clusters()
    cluster_count = len(clusters)
    cluster_map = {clusters[i]: i for i in range(cluster_count) }
    del clusters

    items = [item for item in partition_set.get_all_elements().keys()]
    element_count = len(items)
    item_map = {items[i]: i for i in range(element_count)}
    del items

    # Data handling!

    def __make_sparse_matrix_for_cluster__(cluster: Cluster, _element_count: int, _cluster_count: int) -> lil_matrix:

        result = lil_matrix(np.zeros(shape=(_element_count, _cluster_count), dtype=np.short), shape=(_element_count, _cluster_count))
        #print(_element_count, _cluster_count)
        for item in cluster:
            result[item_map[item], cluster_map[cluster]] = 1
        return result

    def __make_matrix_label_for_cluster__(cluster: Cluster):
        result = np.zeros((len(item_map), len(cluster_map)), dtype=np.short)
        for item in cluster:
            result[item_map[item], cluster_map[cluster]] = 1
        return result.flatten('C').astype(np.short)

    X = []
    y = []

    for partition in tqdm(partition_set):
        for cluster in partition.values():
            label = __make_matrix_label_for_cluster__(cluster)

            for item in cluster:
                X.append(item.as_composition_list())
                y.append(label)

    print("Started ensembling!")

    



    # # Ensembling!
    # print("LR ...")
    # clf1 = LogisticRegression(random_state=1)
    # print("RFC ...")
    # clf2 = RandomForestClassifier(n_estimators=390, random_state=1)
    # print("GNB ..")
    # clf3 = GaussianNB()

    # print("VC Starting...")
    # eclf = VotingClassifier(estimators=[('lr', clf1), ('rfc', clf2), ('gnb', clf3)], voting='hard')
    
    # print("Evaluating! ...")
    # for clf, label in zip([clf1, clf2, clf3, eclf], ['Logistic Regression', 'Random Forest', 'naive Bayes', 'Ensemble']):
    #     scores = cross_val_score(clf, X, y, scoring='accuracy', cv=5, error_score='raise')
    #     print("Accuracy: %0.2f (+/- %0.2f) [%s]" % (scores.mean(), scores.std(), label))


if __name__ == '__main__':
    args = parser()
    run(args)