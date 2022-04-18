from Cluster import Cluster, Partition
from typing import List, Dict, Tuple
from math import factorial, floor, log, sqrt

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

class MRIEvaluator:
    @staticmethod
    def evaluate(eval_partition: Partition, true_partition: Partition, total_object_amount: int) -> float:
        counter, divisor1, divisor2 = 0, 0, 0

        first = True
        for eval_cluster in eval_partition.values():

            divisor1 += Evaluator.__MRI_calc__(len(eval_cluster), total_object_amount)

            for true_cluster in true_partition.values():

                if first:
                    divisor2 += Evaluator.__MRI_calc__(len(true_cluster), total_object_amount)

                intersection_len = len(eval_cluster.intersection(true_cluster))
                if intersection_len <= 0: continue
                counter += Evaluator.__MRI_calc__(total_object_amount * intersection_len, \
                     len(eval_cluster) * len(true_cluster), intersection_len)

            first = False

        divisor = sqrt(divisor1 * divisor2)
        result = counter / divisor
        return result