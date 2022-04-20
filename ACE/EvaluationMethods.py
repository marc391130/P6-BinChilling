from Cluster import Cluster, Partition
from typing import List, Dict, Tuple
from math import factorial, floor, log, sqrt
import PartitionSetReader
import argparse
import os

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
        for eval_cluster in eval_partition.values():
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ACE',
        description="""ACE BINNING ENSEMBLER""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s WRITE THIS LATER PLZ",
        add_help=True
        )
    
    p_args = parser.add_argument_group(title='Setup options', description=None)
    p_args.add_argument('--method', metavar='', required=True, \
        dest='method', help='Methods:\n1: ARI.\n2: NMI')
    p_args.add_argument('--p1', metavar='', required=True, \
        dest='p1', help='Path to result Partition.')
    p_args.add_argument('--p2', metavar='', required=True, \
        dest='p2', help='Path to true Partition.')
    p_args.add_argument('-e', metavar='', required=True, \
        dest='e', help='Number of Objects (Not clusters) in partition [Default = None]')

    args = parser.parse_args()

    p1_path = os.path.abspath(args.p1)
    if os.path.isfile(p1_path) is False:
        raise FileNotFoundError(p1_path)

    p2_path = os.path.abspath(args.p2)
    if os.path.isfile(p2_path) is False:
        raise FileNotFoundError(p2_path)

    number_of_elements = int(args.e)

    p1 = PartitionSetReader.PartitionSetReader.__read_single_partition__(p1_path)
    p2 = PartitionSetReader.PartitionSetReader.__read_single_partition__(p2_path)

    if args.method == '1': print(ARIEvaluator.evaluate(p1, p2, number_of_elements))
    elif args.method == '2': print(NMIEvaluator.evaluate(p1, p2, number_of_elements))