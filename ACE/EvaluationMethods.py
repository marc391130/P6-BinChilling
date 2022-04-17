from Cluster import Cluster, Partition
from typing import List, Dict, Tuple
from math import factorial, floor

class Evaluator:
    def __binomial_coefficient__(self, n) -> float:
        counter = factorial(n)
        divisor = 2 * factorial(n - 2) # Formula is 2! * (n - 2)!, Note 2! is 2.
        result = counter / divisor
        return result


class ARIEvaluator(Evaluator):
    @staticmethod
    def evaluate(self, eval_partition: Partition, true_partition: Partition, total_object_amount: int) -> float:
        nij_comp, ni_comp, nj_comp = 0, 0, 0

        for eval_cluster in eval_partition.values():
            first = True
            ni_comp += self.__binomial_coefficient__(len(eval_cluster))

            for true_cluster in true_partition.values():
                
                cluster_intersection = eval_cluster.intersection(true_cluster)
                nij_comp += self.__binomial_coefficient__(len(cluster_intersection))

                if not first: continue
                nj_comp += self.__binomial_coefficient__(len(true_cluster))
                first = False

        expected_ri = floor(ni_comp*nj_comp) / self.__binomial_coefficient__(total_object_amount)
        counter = nij_comp - expected_ri
        divisor = (0.5 * floor(ni_comp + nj_comp)) - expected_ri

        score = counter / divisor
        return score