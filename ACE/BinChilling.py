from __future__ import annotations
from functools import partialmethod
from multiprocessing import cpu_count
from posixpath import split
from typing import Callable, Dict, List, Tuple, Generic, TypeVar

from AdaptiveEnsembler import Ensembler
from Cluster import Cluster, Partition, PartitionSet
import numpy as np
from tqdm import tqdm
from itertools import islice
import Assertions as Assert
from time import time
from MemberSimularityMatrix import MemberMatrix, MemberSimularityMatrix
from ClusterSimilarityMatrix import SparseClustserSimularity, SparseDictHashMatrix
from AdaptiveEnsemblerExtensions import AssignRegulator, MergeRegulator, QualityMeasuerer, sort_merged_cluster_multithread, sort_merged_cluster_singlethread, __partial_cluster_certainty_degree__, handle_estimate_target_clusters
from AdaptiveEnsemblerDomainExtensions import MergeSCGEvaluator, SCGAssignRegulator
from io import TextIOWrapper
from BinEvaluator import BinEvaluator
from Domain import ContigData

__global_disable_tqdm = False
tqdm.__init__ = partialmethod(tqdm.__init__, disable=__global_disable_tqdm)

THREAD_COUNT = cpu_count()

class MyLogger:
    def __init__(self, console_log: bool = True, logfile: TextIOWrapper = None) -> None:
        self.__should_log__, self.logfile = console_log, logfile 
    
    def __call__(self, string: str) -> None:
        self.log(string)
        
    def log(self, string: str) -> None:
        if self.__should_log__: print(string)
        if self.logfile is not None: print(string, file=self.logfile)

class AbstractEnsembler:
    def __init__(self, logger: MyLogger) -> None:
        self.log = logger
    
    def ensemble(self, gamma: PartitionSet) -> Partition:
        pass
    
    def build_final_partition(self, gamma: PartitionSet, candidate_clusters: List[Cluster]):
        partition = Partition(list(gamma.get_all_elements().keys()))
        found_items, dups = set(), 0
        
        for cluster_index in tqdm(range(len(candidate_clusters))):
            cluster = candidate_clusters[cluster_index]
            for item in cluster:
                if item in found_items:
                    self.log(f'Item {item.name} is a duplicate, having {len(item.SCG_genes)} and len {item.contig_length}')
                    dups += 1
                    continue
                partition.add(str(cluster), item)
                found_items.add(item)
        
        if dups > 0: self.log(f'{dups} duplicate_items')
        self.log(f'{len(found_items)} total items')
        return partition
        
    
class BinChillingEnsembler(AbstractEnsembler):
    def __init__(self, 
                chiller: Chiller,
                binner: Binner,
                bin_eval: BinEvaluator,
                target_clusters_est: int or Callable[[PartitionSet], int] = None,
                chunksize: int = 1,
                processors: int = None,
                logger: MyLogger = None
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
        self.log = logger if logger is not None else MyLogger()
        self.target_clusters_est = target_clusters_est
    
    def ensemble(self, gamma: PartitionSet) -> Partition:
        start_time = time()
        
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
        self.log(f"Finished in time {(time() - start_time):0.02f}s")
        
        return partition
    
    
    def pick_candidate_Clusters(self, cluster_lst: List[Cluster], target_clusters: int,  partition_count: int) -> Tuple[List[Cluster], List[Cluster]]:
        return cluster_lst, []
        decorated_lst = [ (cluster, cluster.mean_member_simularity(partition_count)) for cluster in cluster_lst]
        sort_lst = [x[0] for x in sorted(decorated_lst, key=lambda x: x[1])]
        return sort_lst[:target_clusters], sort_lst[target_clusters:]
    
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
            logger: MyLogger):
        self.a1_min = a1_min
        self.log = logger if logger is not None else MyLogger()
        self.merge_regulator = merge_predicate
        self.alpha_delta = alpha_delta
        self.alpha1 = alpha1
    
    def run(self, gamma: PartitionSet, cluster_lst: List[Cluster], bin_eval: BinEvaluator, processors: int, chunksize: int) -> List[Cluster]:
        i, start_loop, alpha1 = 0, time(), self.alpha1
        item_count = gamma.count_elements()
        
        self.log("Builing cluster similarity matrix...")
        cluster_matrix = SparseClustserSimularity.build_multithread(cluster_lst, item_count, self.a1_min, processors, chunksize)
        
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
                return -1
            elif len(merged_lst) < 10:
                return sort_merged_cluster_singlethread(cluster_matrix, merged_lst)
            return sort_merged_cluster_multithread(cluster_matrix, merged_lst, processors, chunksize)
        
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
                bin_evaluator: BinEvaluator,
                quality_measurer: QualityMeasuerer,
                alpha2: float = 0.75,
                logger: MyLogger = None) -> None:
        self.alpha2 = alpha2
        self.bin_evaluator = bin_evaluator
        self.quality_measure = quality_measurer
        self.log = logger if logger is not None else MyLogger()
    
    
    def assign_item_to_one_cluster(self, gamma: PartitionSet, cluster_lst: List[Cluster],\
        similarity_matrix: MemberSimularityMatrix, non_cand_membermatrix: MemberMatrix) -> Partition:
        
        all_items = gamma.get_all_items()
        
        def sort_by_sim(items: List[Tuple[ContigData, float]]) -> List[ContigData]:
            return [x[0] for x in sorted(items, key=lambda x: x[1], reverse=True )]
        
        def split_lst(items: List[ContigData]) -> Tuple[List[ContigData], List[ContigData]]:
            scg, empty = [], []
            for x in items: (empty if len(x.SCG_genes) == 0 else scg).append(x)
            return scg, empty
        
        self.log("Classifying object certainty...")
        certain_lst, uncertain_lst, lost_lst \
            = self.identify_object_certainy(all_items, similarity_matrix)
        
        self.log("\nFound: ")
        self.log(f"Certain objects: {len(certain_lst)}")
        self.log(f"Uncertain objects: {len(uncertain_lst)}")
        if len(lost_lst) > 0: self.log(f"Lost objects: {len(lost_lst)}")
        self.log('\n')
        
        scg_lst, uncertain_lst = split_lst(sort_by_sim(uncertain_lst))
        
        self.log("Assigning certain objects...")
        candidate_clusters = self.__assign_certains_objects__(certain_lst, cluster_lst)
        
        self.log("Assigning uncertain objects with SCG")        
        candidate_clusters, bad_items = self.__assign_using_SCGs__(scg_lst,\
            similarity_matrix, candidate_clusters, gamma, force=False)

        self.log("Assigning uncertain objects...")
        for item in tqdm(uncertain_lst):
            x =self.__handle_item2__(item, similarity_matrix.get_row(item), similarity_matrix, force=False)
            if x is not None: bad_items.append(x)

        bad_items += lost_lst

        if len(bad_items) > 0:
            self.log(f'Resulted in {len(bad_items)} bad items, trying to relocate using common co-assosication...')
            candidate_clusters = self.assign_remaining(bad_items, candidate_clusters, similarity_matrix, gamma)
        else:
            self.log('No bad items, skipping relocation effert...')

        return candidate_clusters
    
    def assign_remaining(self, item_lst: List[object], cluster_lst: List[Cluster], \
        simularity_matrix: MemberSimularityMatrix, gamma: PartitionSet) -> List[Cluster]:
        
        
        self.log(f'Resulted in {len(item_lst)} bad items, trying to relocate using common co-assosication...')
        self.recalculate_simularity(item_lst, simularity_matrix, gamma, cluster_lst)
        
        sim_lst = [ (x, simularity_matrix.item_max(x)) for x in item_lst ]
        
        def sort_by_sim(items: List[Tuple[ContigData, float]]) -> List[ContigData]:
            return [x[0] for x in sorted(items, key=lambda x: x[1], reverse=True )]
        def split_lst(items: List[ContigData]) -> Tuple[List[ContigData], List[ContigData]]:
            scg, empty = [], []
            for x in items: (empty if len(x.SCG_genes) == 0 else scg).append(x)
            return scg, empty
        
        
        self.log('Re-assigning objects using SCGs...')
        scg_items, rest_lst = split_lst(sort_by_sim(sim_lst))
        
        cluster_lst, bad_items = self.__assign_using_SCGs__(scg_items,\
            simularity_matrix, cluster_lst, gamma, force=True)
        self.log("Re-assigning uncertain objects...")
        for item in tqdm(rest_lst):
            x = self.__handle_item2__(item, simularity_matrix.get_row(item), simularity_matrix, force=True)
            if x is not None: bad_items.append(x)
            
        if len(bad_items) > 0:
            self.log(f'Killing {len(bad_items)} items')
            # self.kill_items(bad_items, cluster_lst)
            # return cluster_lst
            self.log(f'Placing {len(bad_items)} bad items in isolated clusters')
            for bad_item in item_lst:
                isolated_cluster = Cluster()
                isolated_cluster.add(bad_item)
                cluster_lst.append(isolated_cluster)
        
        return cluster_lst
            
        
    
    def recalculate_simularity(self, item_lst: List[object], simularity_matrix: MemberSimularityMatrix, gamma: PartitionSet, cluster_lst: List[Cluster]) \
        -> SparseDictHashMatrix[object, Tuple[float, int]]:
        
        self.log("Building common co-assosiation matrix...")
        all_items = gamma.get_all_items()
        
        member_matrix = MemberMatrix.build(cluster_lst, all_items)
        common_coassosiation = member_matrix.calculate_coassosiation_matrix(set(item_lst), gamma)
        
        for item in tqdm(item_lst):
            for cluster in cluster_lst:
                value = member_matrix.average_common_neighbors(item, cluster, common_coassosiation)
                simularity_matrix[item, cluster] = value
                
        return common_coassosiation
        
    good_clusters = set()
    def __assign_certains_objects__(self, certain_lst: List[Tuple[object, Cluster]], candidate_clusters: List[Cluster]) \
        -> Tuple[List[Cluster], List[Cluster]]:
        
        for item_cluster in tqdm(certain_lst):
            #done so type hinting can actually be done. Damn you tqdm
            item: object = item_cluster[0]
            cluster: Cluster = item_cluster[1]
            
            #remove item from all other clusters in candidate clusters
            for can_cluster in candidate_clusters:
                can_cluster.remove(item)
            #add it back into best cluster
            cluster.add(item)
            self.good_clusters.add(cluster)
        print(f'Found {len(self.good_clusters)} good clusters')
        return candidate_clusters
    
    def kill_items(self, kill_lst: List[ContigData], cluster_lst: List[Cluster]):
        self.log(f'Killing {len(kill_lst)} items')
        for item in tqdm(kill_lst):
            for cluster in cluster_lst:
                cluster.remove(item)
        return cluster_lst
        
    def identify_object_certainy(self, items: List[object], similarity_matrix: MemberSimularityMatrix) \
        -> Tuple[List[Tuple[object, Cluster]], List[Tuple[object, float]], List[object]]: 
            #Totally certain objects,  object with certainty, Totally uncertain objects
        certain_lst, uncertain_lst, lost_lst  = [], [], []
    
        
        for item in tqdm(items):
            cluster, similarity = similarity_matrix.item_argMax(item)
            
            if similarity >= 1:
                certain_lst.append( (item, cluster) )
            elif similarity == 0:
                lost_lst.append( item )
            elif 0 < similarity < 1:
                uncertain_lst.append( (item, similarity) ) 
            else:
                raise Exception(f"something went wrong, simularity value of '{similarity}' outside bound [0, 1]")
        
        return certain_lst, uncertain_lst, lost_lst
        
    
    
    #force indicates whether items should be forcefully placed within a bin or to add it to bad_items return value.
    def __assign_using_SCGs__(self, item_lst: List[ContigData], similarity_matrix: MemberSimularityMatrix, \
        candidate_clusters: List[Cluster],  gamma: PartitionSet, force: bool = False) \
            -> Tuple[List[Cluster], List[ContigData]]: 
        #returns List of certain_clusters, List of negative contigs, 
        bad_items = []
        count = 0
        for item in tqdm(item_lst):
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
                values = self.bin_evaluator.calculate_item_score(cluster, extra_item=item)
                cluster_sim = similarity_matrix.get_column(cluster)
                score1 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items()])
                #score2 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items() if x is not item])
                score2 = score1 - (values[item] * cluster_sim.get(item, 0.0))

                score = similarity * abs(score1 - score2)
                better_with = score1 >= score2
                
                score = score if better_with else score*(-1)

                # score = (similarity) * (score1 - score2)
                print(score)

                if score > best_score:
                    best_score = score
                    best_cluster = cluster
                    
            print("next")
            for cand_cluster in candidate_clusters:
                if cand_cluster is best_cluster: continue #dont remove from best cluster
                if item not in cand_cluster: continue 
                cand_cluster.remove(item)
            
            #add to bad items if none
            if best_cluster is None:
                bad_items.append(item)
                continue    
            
            #ensure item is in cluster
            if item not in best_cluster:
                best_cluster.add(item)
            
            if best_score < 0:
                print('I GOT HERE NIGNOG', force)
                if force:
                    print('Fuck this 2')
                    count += 1
                else:
                    best_cluster.remove(item)
                    bad_items.append(item)
            else: self.good_clusters.add(best_cluster)
            similarity_matrix.assign_item_to(best_cluster, item)
            
        if force: self.log(f'Forcefully placed {count} items in bins')
        print('>REMOVE LATER scgs:' , len(bad_items), ' | good clusters ', len(self.good_clusters))
        return candidate_clusters, bad_items
    
    
    def __handle_item2__(self, item: ContigData, related_clusters: Dict[Cluster, float],\
        similarity_matrix: MemberSimularityMatrix, force: bool = True) -> ContigData or None:
        
        if len(related_clusters) == 0:
            return item #bad
        best_cluster, max_sim = max(related_clusters.items(), key=lambda x: x[1])
        if max_sim > 0.5 or force:
            similarity_matrix.assign_item_to(best_cluster, item)

            if item not in best_cluster:
                best_cluster.add(item)
            
            for cluster, sim in related_clusters.items():
                if cluster is best_cluster: continue
                
                cluster.remove(item)
            return None #i.e., its good
        else:
            return item #i.e., its bad
        
    
    def __handle_item_without_SCGs__(self, item: ContigData, related_clusters: Dict[Cluster, float],\
        gamma: PartitionSet, similarity_matrix: MemberSimularityMatrix\
        ) -> ContigData or None: #Bad items
        
        best_value, best_clusters = np.NINF, []
        for cluster, similarity in related_clusters.items(): #Best clusters
            if similarity > best_value:
                best_clusters = [cluster]
            elif similarity == best_value:
                best_clusters.append(cluster)
        
        if best_value < self.alpha2:
            return item
        
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
        return None
    

