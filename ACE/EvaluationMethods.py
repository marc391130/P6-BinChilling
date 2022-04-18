from Cluster import Cluster, Partition
from typing import List, Dict, Tuple
from math import factorial, floor

class Evaluator:
    @staticmethod
    def __binomial_coefficient__(n) -> float:
        counter = factorial(n)
        divisor = 2 * factorial(n - 2) # Formula is 2! * (n - 2)!, Note 2! is 2.
        result = counter / divisor
        return result


class ARIEvaluator:
    @staticmethod
    def evaluate(eval_partition: Partition, true_partition: Partition, total_object_amount: int) -> float:
        nij_comp, ni_comp, nj_comp = 0, 0, 0

        first = True
        for eval_cluster in eval_partition.values():
            if len(eval_cluster) < 2: continue
            ni_comp += Evaluator.__binomial_coefficient__(len(eval_cluster))

            for true_cluster in true_partition.values():
                if not first or len(true_cluster) < 2: continue
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