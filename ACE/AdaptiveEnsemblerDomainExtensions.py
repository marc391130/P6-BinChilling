from typing import List, Dict, Tuple, Callable
from AdaptiveEnsemblerExtensions import MergeRegulator
from ClusterSimilarityMatrix import SparseClustserSimularity
from Cluster import Cluster, Partition, PartitionSet
from Domain import ContigData
from BinEvaluator import BinEvaluator 
from time import time

class MergeSCGEvaluator(MergeRegulator):
    def __init__(self, a1_min: float, all_SCGs : set, debug: bool = False) -> None:
        self.a1_min = a1_min
        self.bin_evaluator = BinEvaluator(all_SCGs)
        self.debug = debug
        self.log_container = []
        self.LastScore = 0
        self.log_filename = str(time())
        
    
    def set_context(self, gamma: PartitionSet):
        self.__context__ = gamma 
        
    def evaluate(self, alpha1: float, cluster_matrix: SparseClustserSimularity,\
        merged_clusters: List[Cluster[ContigData]]) -> bool:
        if alpha1 < self.a1_min: return self.__log_result__(True, {}, -1, alpha1)

        all_clusters : List[Cluster[ContigData]] = cluster_matrix.get_available_clusters() + merged_clusters
        
        total_dct = self.bin_evaluator.score(all_clusters)
        penalty = len([x for x in total_dct.values() if x <= 0] ) / len(total_dct)
        
        result_value =  (sum(total_dct.values()) / len(total_dct)) * (1 - penalty)
        return_val = self.LastScore > result_value
        
        self.LastScore = result_value
        
        
        # return return_val
        return self.__log_result__(False, total_dct, result_value , alpha1)
        
    
    def __log_result__(self, result: bool, values: Dict[Cluster, float], result_value: float, a1: float) -> bool:
        if not self.debug: result

        if result is True:
            with open('./bin_merge_regulator.txt', 'w') as f:
                for i in range(len(self.log_container)):
                    result_value, zero_values, cluster_count, a1, tup, com_tup = self.log_container[i]
                    near, substantial, moderate, partial, bad = tup
                    completeness, contamination = com_tup
                    f.write( f'{i}> result: {result_value}, zero_values: {zero_values}, total_cluster: {cluster_count} | Comp: {completeness}, cont: {contamination}, avg: { (completeness - contamination) / cluster_count }\n' )
            
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
        
        
        
        
        
        