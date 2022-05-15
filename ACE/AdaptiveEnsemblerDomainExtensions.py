import queue
from sys import maxsize
from typing import List, Dict, Tuple, Callable, Iterator
from EnsemblerTools import MergeRegulator, target_bin_3_4th_count_estimator
from Cluster_matrices import MemberSimularityMatrix, ClustserSimularityMatrix
from ClusterDomain import Cluster, Partition, PartitionSet
from Domain import ContigData
from BinEvaluator import BinEvaluator 
from time import time
import numpy as np

from matplotlib import pyplot as plot

BUFFER_COUNT = 100

class ACEqueue:
    def __init__(self, maxsize) -> None:
        self.maxsize = maxsize
        self.container = []
        self.best_value, self.best_result = np.NINF, []
        
    def put(self, item: Tuple[float, List[Cluster]]):
        if len(self.container) == 0 or item[0] > self.best_value:
            self.best_value, self.best_result = item[0], item[1]
        r = np.NINF
        if self.full():
            r = self.pop()
        self.container.append(item[0])
        return r

    def pop(self) -> float:
        return self.container.pop(0)

    def full(self) -> bool:
        return len(self.container) > self.maxsize

    def get_best_list(self) -> List[Cluster]:
        return self.best_result

    def __iter__(self) -> Iterator:
        return self.container.__iter__()
    
    def __len__(self) -> int:
        return len(self.container)

class MergeSCGEvaluator(MergeRegulator):
    def __init__(self, a1_min: float, bin_evaluator: BinEvaluator, debug: bool = False) -> None:
        super().__init__(a1_min)
        # self.a1_min = a1_min
        self.bin_evaluator = bin_evaluator
        self.debug = debug
        self.log_container, self.log_filename = [], str(time())
        self.buffer = ACEqueue(maxsize=BUFFER_COUNT)
        self.LastScore, self.merge_count = 0, 0
        
    def evaluate(self, alpha1: float, cluster_matrix: ClustserSimularityMatrix, merged_clusters: List[Cluster]) \
        -> Tuple[bool, List[Cluster]]:
        self.merge_count += 1
        partitions_count = len(self.__context__)

        all_clusters: List[Cluster[ContigData]] = cluster_matrix.get_available_clusters() + merged_clusters
        non_merged_clusters = [cluster for cluster in all_clusters if cluster.__partition_id__ is not None]
        all_merged_clusters = [cluster for cluster in all_clusters if cluster.__partition_id__ is None]
        
        total_dct = self.bin_evaluator.score_lst(all_clusters)
        pentalty2 = 100 / len(all_clusters ) #- ( (len(non_merged_clusters) / len(all_clusters)) **2 )
        
        score = sum([value for value in total_dct.values()])
        print(score)
        

        result = super().evaluate(alpha1, cluster_matrix, merged_clusters)
        self.__log_result__(result, total_dct, score, alpha1) 
        return result
        
    def __log_result__(self, result: Tuple[bool, List[Cluster]], values: Dict[Cluster, float], result_value: float, a1: float) -> Tuple[bool, List[Cluster]]:
        if not self.debug: return result
        complete, result_lst = result
        if complete == True:
            with open('./bin_merge_regulator.txt', 'w') as f:
                for i in range(len(self.log_container)):
                    result_value, average_sim, cluster_count, a1, com_tup= self.log_container[i]
                    completeness, contamination, mp = com_tup
                    f.write( f'{i}> result: {result_value}, average_sim: {average_sim}, total_cluster: {cluster_count} | a1: { a1 }, ratio: { cluster_count / a1 if a1 > 0.0 else -1 } | Comp: {completeness}, cont: {contamination}, megabin: { mp }\n' )

            z = self.log_container
            xAxis = [i for i in range(len(self.log_container))]
            plot.plot([x[0] for x in self.log_container], label='score')
            # plot.plot([x for x in score2], label='score2')
            plot.plot([x[4][0] / x[2] for x in self.log_container], label='completeness')
            plot.plot([x[4][1] / x[2] for x in self.log_container], label='contaminatin')
            plot.plot([(x[4][2]) / x[2] for x in self.log_container], label='megabin')
            
            # plot.plot([ 2*((x[4][2]*x[4][0]) / (x[4][2] + x[4][0]) ) for x in self.log_container], label='F1-score')
            # plot.plot([ 3*((x[4][2]*x[4][0]*x[1]) / (x[4][2] + x[4][0] + x[1]) ) for x in self.log_container], label='F1-score * average sim')
            plot.plot([ x[1]*(x[3]**2) for x in self.log_container], label='average_sim * a1^2')
            
            
            plot.plot([(x[4][0] - x[4][1]) / x[2] for x in self.log_container], label='total')
            plot.plot([(x[4][0] / x[4][1]) if x[4][1] > 0 else 0 for x in self.log_container], label='ratio')

            plot.legend()
            plot.show()

            self.log_container = []
        else:
            total_completeness, total_contamination, total_mp = 0, 0, 0
            for cluster in values.keys():
                completeness, contamination, mp = self.bin_evaluator.evaluate(cluster)
                total_completeness += completeness
                total_contamination += contamination
                total_mp += mp
            average_sim = sum([x.mean_member_simularity(25) for x in values.keys()]) / len(values)
            # non_zero_count = len([x for x in values.values() if x > 0]) 
            self.log_container.append( (result_value, average_sim, len(values), a1, (total_completeness, total_contamination, total_mp ) ) )
            
        return result