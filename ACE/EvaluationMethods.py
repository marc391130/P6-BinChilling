from Cluster import Cluster, Partition
from typing import List, Dict, Tuple
from math import factorial, floor, log, sqrt
from Domain import ContigData, bin_size
from PartitionSetReader import PartitionSetReader
from ContigReader import ContigReader
import argparse
import os
from os.path import join
from tqdm import tqdm

class Evaluator:
    @staticmethod
    def __binomial_coefficient__(n) -> float:
        counter = factorial(n)
        divisor = 2 * factorial(n - 2) # Formula is 2! * (n - 2)!, Note 2! is 2.
        result = counter / divisor
        return result

    @staticmethod
    def __MRI_calc__(n: int, all_n: int, multiplier: int = None) -> float:
        if multiplier is None: multiplier = n
        return (multiplier * log(n / all_n))


class ARIEvaluator:
    @staticmethod
    def evaluate(eval_partition: Partition, true_partition: Partition, total_object_amount: int) -> float:
        nij_comp, ni_comp, nj_comp = 0, 0, 0

        # print(len(eval_partition), len(true_partition), total_object_amount)

        first = True
        for eval_cluster in tqdm(eval_partition.values()):
            if len(eval_cluster) < 2: continue
            ni_comp += Evaluator.__binomial_coefficient__(len(eval_cluster))

            for true_cluster in true_partition.values():
                if len(true_cluster) < 2: continue

                if first:
                    nj_comp += Evaluator.__binomial_coefficient__(len(true_cluster))

                cluster_intersection_len = len(eval_cluster.intersection(true_cluster))
                if cluster_intersection_len < 2: continue
                nij_comp += Evaluator.__binomial_coefficient__(cluster_intersection_len)

            first = False

        expected_ri = floor(ni_comp * nj_comp) / Evaluator.__binomial_coefficient__(total_object_amount)
        counter = nij_comp - expected_ri
        divisor = (0.5 * floor(ni_comp + nj_comp)) - expected_ri

        score = counter / divisor
        return score

class NMIEvaluator:
    @staticmethod
    def evaluate(eval_partition: Partition, true_partition: Partition, total_object_amount: int) -> float:
        counter, divisor1, divisor2 = 0, 0, 0

        # print(len(eval_partition), len(true_partition), total_object_amount)

        first = True
        for eval_cluster in eval_partition.values():
            if len(eval_cluster) < 1: continue
            divisor1 += Evaluator.__MRI_calc__(len(eval_cluster), total_object_amount)

            for true_cluster in true_partition.values():

                if first:
                    if len(true_cluster) < 1: continue
                    divisor2 += Evaluator.__MRI_calc__(len(true_cluster), total_object_amount)

                intersection_len = len(eval_cluster.intersection(true_cluster))
                if intersection_len <= 0: continue
                counter += Evaluator.__MRI_calc__(total_object_amount * intersection_len, \
                     len(eval_cluster) * len(true_cluster), intersection_len)

            first = False

        divisor = sqrt(divisor1 * divisor2)
        result = counter / divisor
        return result


def get_all_compare_lst(partitions: List[str]) -> List[Tuple[str, str]]:
    return [(partitions[i], partitions[j]) for i in range(len(partitions)) for j in range(i+1, len(partitions))  ]

def get_double_compare_lst(lst1, lst2) -> List[str]:
    return [ (lst1[i], lst2[j]) for i in range(len(lst1)) for j in range(len(lst2)) if i != j]

def calc_bin_size(cluster: Cluster[ContigData]) -> int:
    return sum([x.contig_length for x in cluster])

def filter_bins(partition: Partition, minsize: int, contig_dct: Tuple[str, ContigData]) -> Partition:
    #return partition
    remove_lst = []

    result = Partition()

    for key, cluster in partition.items():
        for item in cluster:
            result.add(key, contig_dct[item])
    

    for cluster in result.values():
        if bin_size(cluster) < minsize:
            remove_lst.append(cluster)

    for cluster in remove_lst:
        result.remove(cluster)

    return result
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ACE',
        description="""ACE BINNING ENSEMBLER""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s WRITE THIS LATER PLZ",
        add_help=True
        )
    
    p_args = parser.add_argument_group(title='Setup options', description=None)
    p_args.add_argument('--method', metavar='', required=False, \
        dest='method', default='BOTH', choices=('ARI', 'NMI', 'BOTH'), help='Choices: ARI, NMI, BOTH)')
    # p_args.add_argument('--p1', metavar='', required=True, \
    #     dest='p1', help='Path to result Partition.')
    # p_args.add_argument('--p2', metavar='', required=True, \
    #     dest='p2', help='Path to true Partition.')
    p_args.add_argument('-F', metavar='', required=False, default=None, \
        dest='folder', help='The partitions folder to compare all tsv files from' )
    p_args.add_argument('-P', metavar='', required=False, nargs='+', default=None, \
        dest='path', help='The partitions to compare' )
    p_args.add_argument('-e', metavar='', required=True, \
        dest='e', help='Number of Objects (Not clusters) in partition [Default = None]')
    p_args.add_argument('-C', help='Combinations settings', choices=('All', 'Half'), required=False, type=str,\
        default='All', metavar='', dest='comp')
    p_args.add_argument('-M', help='Minimum bin size to include', type=int, required=False,\
        default=0, metavar='', dest='minsize')
    p_args.add_argument('--fasta', help='path to fasta file', type=str, required=True,\
        metavar='', dest='fastafile')
    p_args.add_argument('--cache', help='path to abundance file', type=str, required=True,\
        metavar='', dest='cachefile')

    args = parser.parse_args()
    all_comb = args.comp == 'All'
    if args.path is None and args.folder is None:
        raise argparse.ArgumentError(args.folder, 'must provide either -F or -P')

    partition_paths1 = [join(args.folder, x) for x in os.listdir(args.folder) if x.endswith('.tsv')] if args.folder is not None else []
    partition_paths2 = args.path if args.path is not None else []

    if len(partition_paths1) + len(partition_paths2) <= 1:
        raise argparse.ArgumentError(args.path, 'must have at least 2 partitions')

    number_of_elements = int(args.e)
    path_lst = partition_paths1 + partition_paths2
    contig_dct = ContigReader(args.fastafile, numpy_file=args.cachefile).read_file_fast(args.cachefile)
    partitionDct = {path: filter_bins(PartitionSetReader.__read_single_partition__(path), args.minsize, contig_dct) for path in path_lst}
    tuple_lst = get_all_compare_lst(path_lst) if all_comb else get_double_compare_lst(partition_paths1, partition_paths2)
    ARI_results, NMI_results = [], []



    for p1, p2 in tuple_lst:
        part1, path1 = partitionDct[p1], p1
        part2, path2 = partitionDct[p2], p2
        print('Evaluationg: ', path1, path2)
        if args.method == 'ARI' or args.method == 'BOTH': 
            ari_result = ARIEvaluator.evaluate(part1, part2, number_of_elements)
            ARI_results.append(ari_result), print('>ARI: ', ari_result)
            
        if args.method == 'NMI' or args.method == 'BOTH': 
            nmi_result = NMIEvaluator.evaluate(part1, part2, number_of_elements)
            NMI_results.append(nmi_result), print('>NMI: ' ,nmi_result)
            
            
    print('\n____________________\n')
    if len(ARI_results) > 1:
        print('>ARI: ', 'gns:',sum(ARI_results) / len(ARI_results), 'max:', max(ARI_results) )
        
    if len(NMI_results) > 1:
        print('>NMI: ', 'gns:',sum(NMI_results) / len(NMI_results), 'max:', max(NMI_results) )