import queue
from sys import maxsize
from typing import List, Dict, Tuple, Callable, Iterator
from AdaptiveEnsemblerExtensions import MergeRegulator, target_bin_3_4th_count_estimator
from ClusterSimilarityMatrix import SparseClustserSimularity
from Cluster import Cluster, Partition, PartitionSet
from Domain import ContigData
from BinEvaluator import BinEvaluator 
from time import time

from matplotlib import pyplot as plot


class ACEqueue:
    def __init__(self, maxsize) -> None:
        self.maxsize = maxsize
        self.container = []
        self.best_result: List[Cluster] = []
        
    def put(self, item: Tuple[float, List[Cluster]]):
        if len(self.container) == 0 or item[0] > max(self.container):
            self.best_result = item[1]
        r = -1
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
    def __init__(self, a1_min: float, all_SCGs : set, debug: bool = False) -> None:
        self.a1_min = a1_min
        self.bin_evaluator = BinEvaluator(all_SCGs)
        self.debug = debug
        self.log_container, self.log_filename = [], str(time())
        self.buffer = ACEqueue(maxsize=5)
        self.LastScore, self.merge_count = 0, 0
        
    
    def set_context(self, gamma: PartitionSet):
        self.__context__ = gamma
        return target_bin_3_4th_count_estimator(gamma)
        
    def evaluate(self, alpha1: float, cluster_matrix: SparseClustserSimularity,\
        merged_clusters: List[Cluster[ContigData]]) -> bool:
        self.merge_count += 1
        if alpha1 < self.a1_min: return self.__log_result__(True, {}, -1, 1)

        all_clusters: List[Cluster[ContigData]] = cluster_matrix.get_available_clusters() + merged_clusters
        total_dct = self.bin_evaluator.score(all_clusters)
        zero_count = len([x for x in total_dct.values() if x <= 0])
        pentalty =  zero_count / len(total_dct)
        score = ( sum(total_dct.values()) / (1+  zero_count) ) * (1 - pentalty) - (1 / (self.merge_count))
        
        value = self.buffer.put( (score, all_clusters) )
        
        if self.buffer.full() is False:
            return self.__log_result__(False, total_dct, score, 0)
        
        is_better = value > max(self.buffer)
        r_val = 1 if is_better else 0
        # return return_val
        return self.__log_result__(is_better, total_dct, score , r_val)
    
        
    def get_merge_result(self) -> List[Cluster]:
        self.LastScore, self.merge_count = 0, 0
        self.__context__ = None
        result = self.buffer.get_best_list()
        self.buffer = ACEqueue(maxsize=5)
        return result
        
    def __log_result__(self, result: bool, values: Dict[Cluster, float], result_value: float, a1: float) -> bool:
        if not self.debug: return result
        if result == True:
            print('hello world!', result)
            with open('./bin_merge_regulator.txt', 'w') as f:
                for i in range(len(self.log_container)):
                    result_value, zero_values, cluster_count, a1, tup, com_tup = self.log_container[i]
                    near, substantial, moderate, partial, bad = tup
                    completeness, contamination = com_tup
                    average = (completeness - contamination) / cluster_count
                    f.write( f'{i}> result: {result_value}, zero_values: {zero_values}, total_cluster: {cluster_count} | a1: { a1 } | Comp: {completeness}, cont: {contamination}, avg: { average }\n' )

            plot.plot([i for i in range(len(self.log_container))], [x[0] for x in self.log_container])
            plot.show()

            self.log_container = []
        else:
            near, substantial, moderate, partial, bad = 0, 0, 0, 0, 0
            total_completeness, total_contamination = 0, 0
            for cluster in values.keys():
                completeness, contamination = self.bin_evaluator.__calculate_completeness__(cluster), self.bin_evaluator.__calculate_contamination__(cluster) 
                result1 = self.bin_evaluator.__calculate_sight__(completeness, contamination)
                total_completeness += completeness
                total_contamination += contamination
                if result1 == 'near': near += 1
                elif result1 == 'substantial': substantial += 1
                elif result1 == 'moderate': moderate += 1
                elif result1 == 'partial': partial += 1
                elif result1 == 'bad': bad += 1
                else: print("I fucked up whoops")
            non_zero_count = len([x for x in values.values() if x > 0]) 
            self.log_container.append( (result_value, len(values) - non_zero_count, len(values), a1, (near, substantial, moderate, partial, bad), (total_completeness, total_contamination) ) )
            
        return result
        
        
    # def __log_output__(self, values: Dict[Cluster, float], result_value: float, last_value) -> None:
    #     if not self.debug: return
    #     with open(f'./{self.log_filename}_{self.call_count}.txt', 'x') as f:
    #         f.write(f'result_value: {result_value}, last_value: {last_value}, cluster_count: {len(values)}\n\n\n')
            
    #         for cluster, value in values.items():
    #             f.write(f'{cluster}: {value}\n\n')
        
        
        
        
        
        