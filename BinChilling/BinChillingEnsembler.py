from __future__ import annotations
from functools import partialmethod
from multiprocessing import cpu_count
from typing import Callable, Dict, List, Tuple, Generic, TypeVar

from ClusterDomain import Cluster, Partition, PartitionSet
import numpy as np
from tqdm import tqdm
import Assertions as Assert
from time import time
from Cluster_matrices import MemberMatrix, MemberSimularityMatrix, CoAssosiationMatrix, ClustserSimularityMatrix
from EnsemblerTools import AbstractEnsembler, BinLogger, MergeRegulator, sort_merged_cluster_multithread, sort_merged_cluster_tasks, sort_merged_cluster_singlethread, __partial_cluster_certainty_degree__, handle_estimate_target_clusters
from BinEvaluator import BinEvaluator, BinRefiner
from Domain import ContigData
from math import sqrt, ceil, floor

THREAD_COUNT = cpu_count()

class BinChillingEnsembler(AbstractEnsembler):
    def __init__(self, 
                chiller: Chiller,
                binner: Binner,
                bin_eval: BinEvaluator,
                target_clusters_est: int or Callable[[PartitionSet], int] = None,
                chunksize: int = 1,
                processors: int = None,
                logger: BinLogger = None
                ) -> None:
        Assert.assert_not_none(chiller)
        Assert.assert_not_none(binner)
        if processors is None: processors = cpu_count()
        Assert.assert_in_range(processors, 1, cpu_count())
        self.chiller = chiller
        self.binner = binner
        self.bin_eval = bin_eval
        self.chunksize = max(chunksize, 1)
        self.processors = processors 
        self.log = logger if logger is not None else BinLogger()
        self.target_clusters_est = target_clusters_est
    
    def ensemble(self, gamma: PartitionSet) -> Partition:
        partition_count = len(gamma)
        all_clusters = gamma.get_all_clusters()
        all_items = gamma.get_all_items()
        
        if partition_count == 0 or len(all_clusters) == 0 or len(all_items) == 0:
            raise Exception('Enesemble is missing input. Either no partitions, clusters or items.')
        
        self.log(f'Starting ensemblement using {len(all_items)} items distributed among {len(all_clusters)} clusters and {len(gamma)} partitions.')        
        target_clusters = handle_estimate_target_clusters(gamma, self.target_clusters_est)
        
        self.chiller.merge_regulator.set_context(gamma, target_clusters)
        available_clusters = self.chiller.run(gamma, all_clusters, self.bin_eval, self.processors, self.chunksize)
        candidate_clusters, non_cand_clusters = self.pick_candidate_Clusters(available_clusters, target_clusters, partition_count)
        
        self.log(f"Finished merging process with {len(available_clusters)} clusters and lambda of {target_clusters}")

        #Consider building non cand clusters
        candidate_clusters = [x if x.__partition_id__ is None else self.copy_cluster(x) for x in candidate_clusters]
        non_cand_membermatrix = MemberMatrix.build(non_cand_clusters, all_items)
        simularity_matrix = MemberSimularityMatrix.IndependentBuild(candidate_clusters, gamma)

        final_clusters = self.binner.assign_item_to_one_cluster(gamma,\
            candidate_clusters, simularity_matrix, non_cand_membermatrix)        
        
        self.log("Building final partition from candidate clusters...")
        partition = self.build_final_partition(gamma, final_clusters)
        self.log(f"Found {len(partition)} total clusters")
        
        return partition
    
    
    def pick_candidate_Clusters(self, cluster_lst: List[Cluster], target_clusters: int,  partition_count: int) -> Tuple[List[Cluster], List[Cluster]]:
        return cluster_lst, []
    
    def copy_cluster(self, cluster: Cluster) -> Cluster:
        new_cluster = Cluster()
        for item in cluster:
            new_cluster.add(item)
        return new_cluster

class Chiller:
    def __init__(self, 
            a1_min,
            alpha1: float, 
            merge_predicate: MergeRegulator,
            alpha_delta: float,
            logger: BinLogger):
        self.a1_min = a1_min
        self.log = logger if logger is not None else BinLogger()
        self.merge_regulator = merge_predicate
        self.alpha_delta = alpha_delta
        self.alpha1 = alpha1
    
    def run(self, gamma: PartitionSet, cluster_lst: List[Cluster], bin_eval: BinEvaluator, processors: int, chunksize: int) -> List[Cluster]:
        i, start_loop, alpha1 = 0, time(), self.alpha1
        item_count = gamma.count_elements()
        
        self.log("Builing cluster similarity matrix...")
        cluster_matrix = ClustserSimularityMatrix.build_multithread(cluster_lst, item_count, self.a1_min, processors, chunksize)
        
        def log_loop() -> None:
            nonlocal i, start_loop, alpha1
            self.log(f"iteration {i:05d}, alpha1: {alpha1:.4f} / { self.merge_regulator.a1_min }, clusters: { len(cluster_matrix.__all_clusters__) }, Time: {(time() - start_loop):0.02f}s")
        
        def find_merge_clusters() -> Tuple[List[Cluster], float]:
            merged_lst, skip_set, max_simularity = [], set(), -1 #only a set, so we dont have to convert it later 
            for clusterTuple, similarity in cluster_matrix.get_entries():
                cluster1, cluster2 = clusterTuple
                if cluster1 in skip_set or cluster2 in skip_set: continue
                if similarity >= alpha1:
                    new_cluster = Cluster.merge(cluster1, cluster2)
                    skip_set.add(cluster1)
                    skip_set.add(cluster2)
                    merged_lst.append(new_cluster)
                else: max_simularity = max(max_simularity, similarity) 
            # for cluster in skip_set: cluster_matrix.remove_cluster(cluster)
            cluster_matrix.remove_set_of_clusters(skip_set)
            return (merged_lst, max_simularity)
        
        def sort_merged_clusters(merged_lst: List[Cluster]) -> float:
            if len(merged_lst) == 0:
                self.log('No clusters to merge')
                return -1.0
            elif len(merged_lst) < 10:
                return sort_merged_cluster_singlethread(cluster_matrix, merged_lst)
            return sort_merged_cluster_tasks(cluster_matrix, merged_lst, processors, chunksize)
        
        try:
            while True:
                i = i+1
                log_loop()
                merged_clusters, max_merge_simularity = find_merge_clusters()
                complete, best_result = self.merge_regulator.evaluate(alpha1, cluster_matrix, merged_clusters)
                if complete:
                    return best_result
                
                max_merged_simularity = sort_merged_clusters(merged_clusters)
                alpha1 = max(max_merge_simularity, max_merged_simularity, alpha1 - self.alpha_delta)
                #alpha1 = min(1, alpha1) #make sure float doesnt end up in 1.00000002
                # alpha1 = alpha1 - self.delta_alpha
        finally:
            self.merge_regulator.get_merge_result()
    
class Binner:
    def __init__(self,
            bin_refiner: BinRefiner,
            bin_evaluator: BinEvaluator,
            alpha2: float = 0.75,
            logger: BinLogger = None) -> None:
        self.bin_refiner = bin_refiner
        self.alpha2 = alpha2
        self.bin_evaluator = bin_evaluator
        self.log = logger if logger is not None else BinLogger()
    
    def sort_by_sim(self, items: List[Tuple[ContigData, float]]) -> List[ContigData]:
        return [x[0] for x in sorted(items, key=lambda x: x[1], reverse=True )]
    
    def assign_item_to_one_cluster(self, gamma: PartitionSet, cluster_lst: List[Cluster],\
        similarity_matrix: MemberSimularityMatrix, non_cand_membermatrix: MemberMatrix) -> Partition:
        
        partition_count = len(gamma)
        
        all_items = gamma.get_all_items()
        
        self.log("\nBuilding co-assosiation matrix...")
        co_matrix = CoAssosiationMatrix.build(gamma)

        recalc_lst = all_items
        old_len = len(all_items)
        loop_min_assign = ceil(sqrt(len(all_items)))
        while True:
            self.log('\n')
            bad_scgs, bad_items = [], []
            certain_lst, scg_lst, uncertain_lst, lost_lst =\
            self.identify_object_certainy(recalc_lst, similarity_matrix)

            if len(certain_lst) > 0:
                cluster_lst = self.__assign_certains_objects__(certain_lst, cluster_lst)
                        
            if len(scg_lst) > 0:
                cluster_lst, bad_scgs = self.__assign_using_SCGs__(self.sort_by_sim(scg_lst), similarity_matrix,\
                    cluster_lst, force=False)
            
            if len(uncertain_lst) > 0:
                cluster_lst, bad_items = self.__assign_uncertain_items_noSCG__(self.sort_by_sim(uncertain_lst), cluster_lst,\
                    similarity_matrix, force=False)

            recalc_lst: List[ContigData] = lost_lst + bad_scgs + bad_items
            
            #Estimate how many to isolate
            assignment_count = old_len - len(recalc_lst)
            isolate_count = loop_min_assign - assignment_count

            if isolate_count > 0 and len(recalc_lst) != 0:
                sorted_recalc = sorted(recalc_lst, key=lambda x: x.contig_length, reverse=True)
                isolate_lst, recalc_lst = sorted_recalc[:isolate_count], sorted_recalc[isolate_count:]
                cluster_lst = self.isolate_items(isolate_lst, cluster_lst)
            if len(recalc_lst) == 0:
                break
            else:
                self.remove_empty_clusters(cluster_lst, similarity_matrix)
                self.recalculate_simularity_fast(recalc_lst, similarity_matrix, cluster_lst, co_matrix, partition_count)
                self.log(f"Managed to assign {old_len - len(recalc_lst)} items, {len(recalc_lst)} items remain...")
                old_len = len(recalc_lst)
                
            
        #loop break
        cluster_lst = self.remove_empty_clusters(cluster_lst, similarity_matrix)
        cluster_lst = self.bin_refiner.Refine(cluster_lst, co_matrix)
        
        return cluster_lst
        

    def identify_object_certainy(self, items: List[ContigData], similarity_matrix: MemberSimularityMatrix) \
        -> Tuple[List[Tuple[object, Cluster]], List[Tuple[object, float]], List[Tuple[object, float]], List[object]]: 
            #Totally certain objects,  object with certainty, Totally uncertain objects
        certain_lst, scg_lst, uncertain_lst, lost_lst  = [], [], [], []
    
        
        for item in tqdm(items, desc='Measuring contigs certainty', bar_format=''):
            item: ContigData
            cluster, similarity = similarity_matrix.item_argMax(item)
            
            if similarity >= 1:
                certain_lst.append( (item, cluster) )
            elif similarity == 0:
                lost_lst.append( item )
            elif 0 < similarity < 1:
                if len(item.SCG_genes) > 0: scg_lst.append( (item, similarity) )
                else: uncertain_lst.append( (item, similarity) ) 
            else:
                raise Exception(f"something went wrong, simularity value of '{similarity}' outside bound [0, 1]")
        
        return certain_lst, scg_lst, uncertain_lst, lost_lst
    
    def __assign_certains_objects__(self, certain_lst: List[Tuple[object, Cluster]], candidate_clusters: List[Cluster]) \
        -> List[Cluster]:
        
        for item_cluster in tqdm(certain_lst, desc='Assigning certain contigs  '):
            #done so type hinting can actually be done. Damn you tqdm
            item: object = item_cluster[0]
            cluster: Cluster = item_cluster[1]
            
            #remove item from all other clusters in candidate clusters
            self.remove_from_all(item, candidate_clusters)
            #add it back into best cluster
            cluster.add(item)
        return candidate_clusters
    
    #force indicates whether items should be forcefully placed within a bin or to add it to bad_items return value.
    def __assign_using_SCGs__(self, item_lst: List[ContigData], similarity_matrix: MemberSimularityMatrix, \
        cluster_lst: List[Cluster], force: bool = False) \
            -> Tuple[List[Cluster], List[ContigData]]: 
        #returns List of certain_clusters, List of negative contigs, 
        bad_items = []
        count = 0
        for item in tqdm(item_lst, desc='Assigning contigs with SCGs'):
            item: ContigData
            #Handle if no SCG
            if len(item.SCG_genes) == 0:
                Assert.assert_fail(f'Item {item} is trying to assign based on SCGs but has none')
            
            row_data = similarity_matrix.get_row(item)
            #if item has no simularity, add to bad and skip to next
            if len(row_data) == 0:
                bad_items.append(item)
                continue
            
            best_cluster: Cluster = None
            best_score: float = np.NINF

            #Try and assign to bin with best impact.
            #Impact being defined as similarity * score_diff
            for cluster, similarity in row_data.items():
                if len(cluster) == 0: continue
                if similarity <= 0: continue
                # values = self.bin_evaluator.score_items(cluster, extra_item=item)
                # cluster_sim = similarity_matrix.get_column(cluster)
                # score1 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items()])
                # #score2 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items() if x is not item])
                # score2 = score1 - (values[item] * cluster_sim.get(item, 0.0))

                score1 = self.bin_evaluator.score(cluster, include_item=item)
                score2 = self.bin_evaluator.score(cluster, skip_item=item)

                score = similarity * (score1 - score2)
                # score = similarity * (score1 if score1 >= score2 else score1 - score2 )

                if score > best_score:
                    best_score = score
                    best_cluster = cluster
                    
            # print("next, best score: ", best_score)
            self.remove_from_all(item, cluster_lst)
            
            #add to bad items if none
            if best_cluster is None:
                bad_items.append(item)
                continue    
            
            #ensure item is in cluster
            if item not in best_cluster:
                best_cluster.add(item)
            
            if best_cluster.__partition_id__ is not None: 
                Assert.assert_fail('A cluster is a source cluster, this shouldnt be possible and will influence the result')
            
            if best_score < 0.0:
                if force:
                    count += 1
                else:
                    best_cluster.remove(item)
                    bad_items.append(item)
            similarity_matrix.assign_item_to(best_cluster, item)
            
        if force: self.log(f'Forcefully placed {count} items in bins')
        return cluster_lst, bad_items
    
    def remove_from_all(self, item: ContigData, cluster_lst: List[Cluster]) -> None:
        for cluster in cluster_lst:
            cluster.remove(item)
    
    def __assign_uncertain_items_noSCG__(self, item_lst: List[ContigData], cluster_lst: List[Cluster],
        similarity_matrix: MemberSimularityMatrix, force: bool = True) -> Tuple[List[Cluster], List[ContigData]]:
        # return self.__assign_using_dna_len__(item_lst, cluster_lst, similarity_matrix)
    
        bad_items = []
        for item in tqdm(item_lst, desc='Assigning uncertain contigs'):
            related_clusters = similarity_matrix.get_row(item)
            if len(related_clusters) == 0:
                bad_items.append(item)
                continue
                
            best_cluster, max_sim = max(related_clusters.items(), key=lambda x: x[1])
            if force:
                bad = False
                for cls, value in related_clusters.items():
                    if cls is best_cluster: continue
                    if value >= max_sim: 
                        bad = True
                        break
                if bad: bad_items.append(item)
                break
                        
            
            if max_sim > 0.5 or force:
                self.remove_from_all(item, cluster_lst)
                best_cluster.add(item)
                similarity_matrix.assign_item_to(best_cluster, item)
            else:
                bad_items.append(item)
        
        return cluster_lst, bad_items
    
    def __assign_using_dna_len__(self, item_lst: List[ContigData], cluster_lst: List[Cluster],
        similarity_matrix: MemberSimularityMatrix) -> Tuple[List[Cluster], List[ContigData]]:

        bad_items = []        
        for item in tqdm(item_lst):
            item: ContigData
            related_clusters = sorted(similarity_matrix.get_row(item).items(), key=lambda x: x[1])
            if len(related_clusters) == 0:
                bad_items.append(item)
                continue
            pass
            
            best_score: float  = np.NINF
            best_cluster: Cluster = None
            best_sim: float = 0.0
            
            for cluster, sim in related_clusters:
                if len(cluster) == 0 or sim < 0.2: continue
                score1 = self.bin_evaluator.score_len(cluster, include_item=item)
                score2 = self.bin_evaluator.score_len(cluster, skip_item=item)
                score = sim * (score1 - score2)
                
                if best_score > score:
                    best_cluster = cluster
                    best_score = score
                    best_sim = sim
            #end loop

            self.remove_from_all(item, cluster_lst)
            similarity_matrix.assign_item_to(cluster, item)
            if best_cluster is None:
                bad_items.append(item)
                continue
            
            if item not in best_cluster:
                best_cluster.append(item)
            
            if best_score < best_sim:
                best_cluster.remove(item)
                bad_items.append(item)
                
        #end loop
        return cluster_lst, item_lst
    
    
    def recalculate_simularity_fast(self, item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
        cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix, partition_count: int):
        member_matrix = MemberMatrix.build(cluster_lst, item_lst)
        
        for item in tqdm(item_lst, desc='Recalc. similarity matrix  '):
            for cluster in member_matrix.get_cluster_set(co_matrix.get_row(item).keys()):
                if len(cluster) == 0: 
                    simularity_matrix[item, cluster] = 0.0
                    continue
                value = min(co_matrix.cluster_mean(item, cluster), 1.0)
                simularity_matrix[item, cluster] = value
        return cluster_lst
    
    def recalculate_simularity(self, item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
        cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix, partition_count: int):
        
        for item in tqdm(item_lst, desc='Recalc. similarity matrix  '):
            for cluster in cluster_lst:
                if len(cluster) == 0: 
                    simularity_matrix[item, cluster] = 0.0
                    continue
                value = min(co_matrix.cluster_mean(item, cluster), 1.0)
                simularity_matrix[item, cluster] = value if 0.0 < value <= 1.0 else 0
        return cluster_lst
    
    def isolate_items(self, item_lst: List[ContigData], cluster_lst: List[Cluster]) -> List[Cluster]:
        for item in item_lst:
            self.remove_from_all(item, cluster_lst)
            cluster = Cluster()
            cluster.add(item)
            cluster_lst.append(cluster)
        return cluster_lst

    def remove_empty_clusters(self, cluster_lst: List[Cluster], similarity_matrix: MemberSimularityMatrix, should_log:bool = False) -> List[Cluster]:
        remove_lst = []
        for cluster in cluster_lst:
            if len(cluster) == 0:
                remove_lst.append(cluster)
                similarity_matrix.remove_cluster(cluster)
        
        if should_log:
            self.log(f'Found {len(remove_lst)} empty clusters...')
        for r_cls in remove_lst:
            cluster_lst.remove(r_cls)
        return cluster_lst