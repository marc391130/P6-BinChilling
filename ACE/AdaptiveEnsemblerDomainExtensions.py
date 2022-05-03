import queue
from sys import maxsize
from typing import List, Dict, Tuple, Callable, Iterator
from AdaptiveEnsemblerExtensions import MergeRegulator, target_bin_3_4th_count_estimator, QualityMeasuerer, AssignRegulator
from MemberSimularityMatrix import MemberSimularityMatrix
from ClusterSimilarityMatrix import SparseClustserSimularity
from Cluster import Cluster, Partition, PartitionSet
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
        self.a1_min = a1_min
        self.bin_evaluator = bin_evaluator
        self.debug = debug
        self.log_container, self.log_filename = [], str(time())
        self.buffer = ACEqueue(maxsize=BUFFER_COUNT)
        self.LastScore, self.merge_count = 0, 0
        
    
    def set_context(self, gamma: PartitionSet, target_clustsers: int) -> None:
        self.__context__ = gamma
        self.target_clusters = target_clustsers
        
    def evaluate(self, alpha1: float, cluster_matrix: SparseClustserSimularity,\
        merged_clusters: List[Cluster[ContigData]]) -> Tuple[bool, List[Cluster]]:
        self.merge_count += 1
        if alpha1 < self.a1_min: return (self.__log_result__(True, {}, -1, alpha1), self.buffer.get_best_list())
        partitions_count = len(self.__context__)

        all_clusters: List[Cluster[ContigData]] = cluster_matrix.get_available_clusters() + merged_clusters
        non_merged_clusters = [cluster for cluster in all_clusters if cluster.__partition_id__ is not None]
        merged_clusters = [cluster for cluster in all_clusters if cluster.__partition_id__ is None]
        
        # xxx_dct = self.bin_evaluator.score(non_merged_clusters)
        total_dct = self.bin_evaluator.score(all_clusters)
        # zero_count = len([x for x in total_dct.values() if x <= 0])
        # pentalty =  zero_count / len(total_dct)
        pentalty2 = 100 / len(all_clusters ) #- ( (len(non_merged_clusters) / len(all_clusters)) **2 )
        
        # value_lst = [(value, cluster.mean_member_simularity(partitions_count)**2) for cluster, value in total_dct.items()]
        #score = sum([value for cluster, value in total_dct.items()])
        score = self.merge_count
        #harmonic mean
        # divisor = sum([(x[1] / x[0] if x[0] != 0 else 0) for x in value_lst])
        # score = sum([x[1] for x in value_lst]) / divisor if divisor != 0 else 0
        # values_lst = [(self.bin_evaluator.__calculate_completeness__(cluster), self.bin_evaluator.__calculate_purity__(cluster), cluster.mean_member_simularity(partitions_count)) for cluster in all_clusters]
        #score = 3 * (sum([x[0] * x[1] * x[2] for x in values_lst]) / sum([x[0] + x[1] + x[2] for x in values_lst]))
        #score *= pentalty2
        print(score)
        
        
        value = self.buffer.put( (score, all_clusters) )
        
        if self.buffer.full() is False:
            return (self.__log_result__(False, total_dct, score, alpha1), self.buffer.get_best_list())
        
        complete = value > max(self.buffer)
        
        # return return_val
        print('>Best value', self.buffer.best_value)
        return (self.__log_result__(complete, total_dct, score , alpha1), self.buffer.get_best_list())
    
        
    def get_merge_result(self) -> List[Cluster]:
        self.LastScore, self.merge_count = 0, 0
        self.__context__ = None
        result = self.buffer.get_best_list()
        self.buffer = ACEqueue(maxsize=BUFFER_COUNT)
        return result
        
    def __log_result__(self, result: bool, values: Dict[Cluster, float], result_value: float, a1: float) -> bool:
        if not self.debug: return result
        if result == True:
            score2 = []
            with open('./bin_merge_regulator.txt', 'w') as f:
                for i in range(len(self.log_container)):
                    result_value, average_sim, cluster_count, a1, com_tup, sco_tup = self.log_container[i]
                    completeness, contamination, purity = com_tup
                    near, substanitial, moderate, partial, bad = sco_tup
                    score2.append( (5*near + 4*substanitial + 3* moderate + 2*partial + bad) /cluster_count  )
                    f.write( f'{i}> result: {result_value}, average_sim: {average_sim}, total_cluster: {cluster_count} | a1: { a1 }, ratio: { cluster_count / a1 if a1 > 0.0 else -1 } | Comp: {completeness}, cont: {contamination}, purity: { purity } | near: {near}, sub: {substanitial}, mode: {moderate}, part: {partial}, bad: {bad}\n' )

            z = self.log_container
            xAxis = [i for i in range(len(self.log_container))]
            plot.plot([x[0] for x in self.log_container], label='score')
            # plot.plot([x for x in score2], label='score2')
            plot.plot([x[4][0] / x[2] for x in self.log_container], label='completeness')
            plot.plot([x[4][1] / x[2] for x in self.log_container], label='contaminatin')
            plot.plot([(x[4][2]*100) / x[2] for x in self.log_container], label='purity')
            
            # plot.plot([ 2*((x[4][2]*x[4][0]) / (x[4][2] + x[4][0]) ) for x in self.log_container], label='F1-score')
            # plot.plot([ 3*((x[4][2]*x[4][0]*x[1]) / (x[4][2] + x[4][0] + x[1]) ) for x in self.log_container], label='F1-score * average sim')
            plot.plot([ x[1]*(x[3]**2) for x in self.log_container], label='average_sim * a1^2')
            
            
            plot.plot([(x[4][0] - x[4][1]) / x[2] for x in self.log_container], label='total')
            plot.plot([(x[4][0] / x[4][1]) if x[4][1] > 0 else 0 for x in self.log_container], label='ratio')
            plot.plot([(x[1] * x[4][0] / x[4][1]) if x[4][1] > 0 else 0 for x in self.log_container], label='avg_sim * ratio')
            # plot.plot([z[i][0]* (self.target_clusters / z[i][2]) for i in range(len(self.log_container)) ], label='Marcus *score')
            # plot.plot([0.1*i* (1 - (self.target_clusters / z[i][2])) for i in range(len(self.log_container)) ], label='Marcus 0.1*score')
            # plot.plot([(1 - (self.target_clusters / z[i][2])) for i in range(len(self.log_container)) ], label='Marcus 1 -(lambda/c) score')
            # plot.plot([x[3] for x in self.log_container], label='alpha1')
            plot.legend()
            #plot.show()

            self.log_container = []
        else:
            total_completeness, total_contamination, total_purity = 0, 0, 0
            near, substanitial, moderate, partial, bad = 0, 0, 0, 0, 0
            for cluster in values.keys():
                scgs = self.bin_evaluator.__calculate_number_of_SCGs__(cluster)
                completeness, contamination, purity = self.bin_evaluator.__calculate_completeness__(scgs), self.bin_evaluator.__calculate_contamination__(scgs) , self.bin_evaluator.__calculate_purity__(scgs)
                result1 = self.bin_evaluator.__calculate_sight__(completeness, contamination)
                if result1 == 'near': near+=1
                elif result1 == 'substanitial': substanitial+=1
                elif result1 == 'moderate': moderate+=1
                elif result1 == 'partial': partial+=1
                elif result1 == 'bad': bad+=1
                
                total_completeness += completeness
                total_contamination += contamination
                total_purity += purity
            average_sim = sum([x.mean_member_simularity(25) for x in values.keys()]) / len(values)
            # non_zero_count = len([x for x in values.values() if x > 0]) 
            self.log_container.append( (result_value, average_sim, len(values), a1, (total_completeness*10000 / len(values), total_contamination*100 / len(values), total_purity*100 / len(values)), (near, substanitial, moderate, partial, bad) ) )
            
        return result
        
        
    # def __log_output__(self, values: Dict[Cluster, float], result_value: float, last_value) -> None:
    #     if not self.debug: return
    #     with open(f'./{self.log_filename}_{self.call_count}.txt', 'x') as f:
    #         f.write(f'result_value: {result_value}, last_value: {last_value}, cluster_count: {len(values)}\n\n\n')
            
    #         for cluster, value in values.items():
    #             f.write(f'{cluster}: {value}\n\n')
class SCGAssignRegulator(AssignRegulator):
    def __init__(self, quality_measure: QualityMeasuerer, merge_regulator: MergeSCGEvaluator, logger: Callable[[str], None] = None) -> None:
        super().__init__(quality_measure, logger)
        self.bin_evaluator = merge_regulator.bin_evaluator

    def assign_items(self, candidate_clusters: List[Cluster], totally_certain_lst: List[Tuple[object, Cluster]], certain_lst: List[Tuple[object, Cluster]], \
        uncertain_lst: List[object], totally_uncertain_map: Dict[object, object], gamma: PartitionSet, similarity_matrix:MemberSimularityMatrix, lost_items: List[object]) -> List[Cluster]:

        self.log("Assigning totally certain objects...")
        candidate_clusters = self.__assign_certains__(totally_certain_lst, candidate_clusters)

        certain_lst = sorted([x for x, y in certain_lst] + uncertain_lst, key=lambda x: similarity_matrix.item_mean(x), reverse=True )

        self.log("Assign Certain Objects...")
        candidate_clusters = self.__assign_using_SCGs__(certain_lst, similarity_matrix, candidate_clusters, gamma)
        # for item in certain_lst:
        #     for cluster in candidate_clusters:
        #         cluster.remove(item)
        
        # self.log("Assign uncertain objects...")
        # candidate_clusters = self.__handle_SCG_certain__(uncertain_lst, similarity_matrix, candidate_clusters, gamma)
        # # for item in uncertain_lst:
        # #     for cluster in candidate_clusters:
        # #         cluster.remove(item)

        self.log("handling lost objects...")
        candidate_clusters = self.__assign_lost_objects__(candidate_clusters, totally_uncertain_map, lost_items)

        return candidate_clusters

    def __assign_using_SCGs__(self, item_lst: List[ContigData], similarity_matrix: MemberSimularityMatrix, candidate_clusters: List[Cluster], \
        gamma: PartitionSet) -> Tuple[List[Cluster], List[ContigData]]: 
        #List of cleaned clusters, List of badly placed contigs
        count = 0
        for item in item_lst:
            row_data = similarity_matrix.get_row(item)
            if len(item.SCG_genes) == 0:
                self.__handle_item_without_SCGs__(item, row_data, gamma, similarity_matrix)
                print('skipping empty contig')
                continue

            best_cluster: Cluster = None
            best_score: float = np.NINF

            for cluster, similarity in row_data.items():
                if len(cluster) == 0: continue
                #score1 = self.bin_evaluator.calculate_score(cluster)
                #score2 = self.bin_evaluator.calculate_score(cluster, item)
                values = self.bin_evaluator.calculate_item_score(cluster, extra_item=item)
                cluster_sim = similarity_matrix.get_column(cluster)
                score1 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items()])
                score2 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items() if x is not item])

                score = (similarity) * (score1 - score2)
                print(score)

                if score > best_score:
                    best_score = score
                    best_cluster = cluster
            print('>next')
            if best_score < 0:
                count += 1 
                #best_cluster.remove(item)
            for cand_cluster in candidate_clusters:
                if cand_cluster is best_cluster: continue
                if item not in cand_cluster: continue 
                cand_cluster.remove(item)

            similarity_matrix.assign_item_to(best_cluster, item)
        print('>Found bad: ', count)
        return candidate_clusters
    
    def __handle_item_without_SCGs__(self, item: ContigData, related_clusters: Dict[Cluster, float], gamma: PartitionSet, similarity_matrix: MemberSimularityMatrix) -> None:
        
        best_value = np.NINF
        best_clusters = []
        for cluster, similarity in related_clusters.items():
            if similarity > best_value:
                best_clusters = [cluster]
            elif similarity == best_value:
                best_clusters.append(cluster)
        
        best_cluster = None
        if len(best_clusters) > 1:
            biggest_change = np.NINF
            for cluster in best_clusters:
                quality_change: float = self.quality_measure.calculate_excluding_quality(cluster, item, len(gamma) - 1, similarity_matrix)
                if quality_change > biggest_change:
                    best_cluster = cluster
                    biggest_change = quality_change
        else:
            best_cluster = best_clusters[0]

        similarity_matrix.assign_item_to(best_cluster, item)
        for cluster in related_clusters.keys():
            if cluster is not best_cluster:
                cluster.remove(item)
