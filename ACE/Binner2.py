import itertools
from os import remove

from math import sqrt, ceil
from BinChilling import MyLogger
from BinEvaluator import BinEvaluator
from BinRefiner import BinRefiner
from Cluster import PartitionSet, Partition, Cluster
from Domain import ContigData
from typing import List, Dict, Tuple, TypeVar, Generic
from MemberSimularityMatrix import CoAssosiationMatrix, MemberMatrix, MemberSimularityMatrix
from tqdm import tqdm
import Assertions as Assert
import numpy as np

from SparseMatrix_implementations import SortKeysByHash, SparseDictHashMatrix

class Binner2:
    def __init__(self,
            bin_refiner: BinRefiner,
            bin_evaluator: BinEvaluator,
            alpha2: float = 0.75,
            logger: MyLogger = None) -> None:
        self.bin_refiner = bin_refiner
        self.alpha2 = alpha2
        self.bin_evaluator = bin_evaluator
        self.log = logger if logger is not None else MyLogger()
    
    def sort_by_sim(self, items: List[Tuple[ContigData, float]]) -> List[ContigData]:
        return [x[0] for x in sorted(items, key=lambda x: x[1], reverse=True )]
    
    def assign_item_to_one_cluster(self, gamma: PartitionSet, cluster_lst: List[Cluster],\
        similarity_matrix: MemberSimularityMatrix, non_cand_membermatrix: MemberMatrix) -> Partition:
        
        # cluster_lst2 = [x for x in cluster_lst if x.__partition_id__ is None]
        # print(f'Yeeted {len(cluster_lst) - len(cluster_lst2)}')
        # cluster_lst = cluster_lst2
        partition_count = len(gamma)
        
        self.log("Classifying object certainty...")
        all_items = gamma.get_all_items()
        
        self.log("Building co-assosiation matrix...")
        co_matrix = CoAssosiationMatrix.build(gamma)

        recalc_lst = all_items
        old_len = len(all_items)
        
        while True:
            self.log('\nIdentify item certainy...')
            certain_lst, scg_lst, uncertain_lst, lost_lst =\
            self.identify_object_certainy(recalc_lst, similarity_matrix)
            
            # self.log('Assigning certain items...')
            cluster_lst = self.__assign_certains_objects__(certain_lst, cluster_lst)
                        
            # self.log('Assigning uncertain items with SCG items...')
            cluster_lst, bad_scgs = self.__assign_using_SCGs__(self.sort_by_sim(scg_lst), similarity_matrix,\
                cluster_lst, force=False)
            
            # self.log('Assigning uncertain items...')
            cluster_lst, bad_items2 = self.__assign_uncertain_items_noSCG__(self.sort_by_sim(uncertain_lst), cluster_lst,\
                similarity_matrix, force=False)

            recalc_lst: List[ContigData] = lost_lst + bad_scgs + bad_items2 #+ [x[0] for x in scg_lst]
            
            #Isolate a part of the clusters if none was placed
            if len(recalc_lst) == old_len and len(recalc_lst) != 0:
                isolate_count = max(ceil(sqrt(len(recalc_lst))), 1)
                sorted_recalc = sorted(recalc_lst, key=lambda x: x.contig_length, reverse=True)
                isolate_lst, recalc_lst = sorted_recalc[:isolate_count], sorted_recalc[isolate_count:]
                cluster_lst = self.isolate_items(isolate_lst, cluster_lst)

            self.log(f"Managed to assign {old_len - len(recalc_lst)} items, {len(recalc_lst)} items remain...")
            if len(recalc_lst) == 0:
                break
            else:
                self.remove_empty_clusters(cluster_lst, similarity_matrix)
                self.recalculate_simularity(recalc_lst, similarity_matrix, cluster_lst, co_matrix, partition_count)
                old_len = len(recalc_lst)
            
        #loop break
        # cluster_lst += self.Handle_remaining_items(recalc_lst, co_matrix)
        cluster_lst = self.remove_empty_clusters(cluster_lst, similarity_matrix)
        # cluster_lst = self.kill_items(recalc_lst, cluster_lst)
        # cluster_lst = self.isolate_items(recalc_lst, cluster_lst)
        # remaining_lst = self.Handle_remaining_items(recalc_lst, gamma)
        # cluster_lst += remaining_lst
        
        #cluster_lst = self.refine_bins(cluster_lst, co_matrix)
        cluster_lst = self.bin_refiner.Refine(cluster_lst, co_matrix)
        
        return cluster_lst
        

    def identify_object_certainy(self, items: List[ContigData], similarity_matrix: MemberSimularityMatrix) \
        -> Tuple[List[Tuple[object, Cluster]], List[Tuple[object, float]], List[Tuple[object, float]], List[object]]: 
            #Totally certain objects,  object with certainty, Totally uncertain objects
        certain_lst, scg_lst, uncertain_lst, lost_lst  = [], [], [], []
    
        
        for item in tqdm(items):
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
        
        for item_cluster in tqdm(certain_lst):
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
                values = self.bin_evaluator.score_items(cluster, extra_item=item)
                cluster_sim = similarity_matrix.get_column(cluster)
                score1 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items()])
                #score2 = sum([y * cluster_sim.get(x, 0.0) for x, y in values.items() if x is not item])
                score2 = score1 - (values[item] * cluster_sim.get(item, 0.0))

                score = similarity * (score1 - score2)
                # score = similarity * (score1 if score1 >= score2 else score1 - score2 )
                # if score1 > score2 and score1 < 0:
                #     score *= -1
                    
                # if score1 >= score2:
                #     score = similarity * score1
                # else:
                #     score = similarity * (score1 - score2)
                    
                #print('sanity score:',score, '| score1: ', score1,'| score2:', score2, '| sim: ', similarity)
                
                # better_with = score1 >= score2
                
                # score = score if better_with else score*(-1)

                # score = (similarity) * (score1 - score2)
                # print(score)

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
        for item in tqdm(item_lst):
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
    
    
    def recalculate_simularity(self, item_lst: List[object], simularity_matrix: MemberSimularityMatrix,\
        cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix, partition_count: int):
        
        self.log("Recalculating simularity matrix")
        for item in tqdm(item_lst):
            for cluster in cluster_lst:
                if len(cluster) == 0: 
                    simularity_matrix[item, cluster] = 0.0
                    continue
                value = min(co_matrix.cluster_mean(item, cluster), 1.0)
                simularity_matrix[item, cluster] = value if 0.0 < value <= 1.0 else 0
        return cluster_lst
        
                
    def kill_items(self, kill_lst: List[ContigData], cluster_lst: List[Cluster]):
        self.log(f'Killing {len(kill_lst)} items')
        # return self.isolate_items(kill_lst, cluster_lst)
        for item in tqdm(kill_lst):
            item: ContigData
            print(f"Killing item '{item.name }' with {len(item.SCG_genes)} scgs and len {item.contig_length}")
            for cluster in cluster_lst:
                cluster.remove(item)
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