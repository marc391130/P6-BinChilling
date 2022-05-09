from typing import Callable, List, Dict, Tuple, TypeVar, Iterator, Iterable
from BinEvaluator import BinEvaluator
from Cluster import Cluster
from Domain import ContigData
from MemberSimularityMatrix import CoAssosiationMatrix
from SparseMatrix_implementations import SparseDictHashMatrix, SortKeysByHash
import numpy as np
from tqdm import tqdm
from BinChilling import MyLogger
import itertools
import Assertions as Assert


score_dct: Dict[Cluster, float] = {}

def split_generator(source: List[Cluster], start: int, stop: int or None = None) -> Iterator[Cluster]:
    stop = stop if stop is not None else len(source)
    if start > stop: Assert.assert_fail(f'start is higher or equal to stop ({start} >= {stop})')
    
    return (source[i] for i in range(start, stop))

class BinRefiner:
    def __init__(self, bin_evaluator: BinEvaluator, min_threshold: float, logger: MyLogger) -> None:
        self.bin_evaluator = bin_evaluator
        self.log = logger
        self.min_threshold = min_threshold
        self.score_dct: Dict[Cluster, float] = {}
    
    def __build_common_co_matrix__(self, cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix,\
        old_matrix: SparseDictHashMatrix[Cluster, float] = None) -> SparseDictHashMatrix[Cluster, float]:
        avg_co_matrix = SparseDictHashMatrix[Cluster, float](SortKeysByHash, default_value=0.0)
        
        for i in tqdm(range(len(cluster_lst))):
            c1 = cluster_lst[i]
            for j in range(i+1, len(cluster_lst)):
                c2 = cluster_lst[j]
                value =  co_matrix.bin_mean(c1, c2) if old_matrix is None or old_matrix.has_entry(c1, c2) is False\
                    else old_matrix[c1, c2]
                avg_co_matrix.set_entry(c1, c2, value if value > self.min_threshold else 0.0)
        return avg_co_matrix
        
    def Refine(self, cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix) -> List[Cluster]:
        
        self.log('Calculating cluster co-assosiation matrix to refine bins...')
        value_matrix = self.__build_common_co_matrix__(cluster_lst, co_matrix)
        
        try:
            self.log('Refining bins...')
            source_clusters = cluster_lst    
            while True:
                skip_set = set()
                new_clusters = self.__refine_clusters__(cluster_lst, value_matrix, skip_set)
                if len(new_clusters) == 0:
                    break 
                else:
                    source_clusters = new_clusters
                    self.log(f'Refined {len(skip_set)} clusters into {len(new_clusters)} new clusters.')
                    self.log('Adjusting cluster co-assosiation matrix...')
                
                for cls in skip_set:
                    cluster_lst.remove(cls)
                for cls in new_clusters:
                    cluster_lst.append(cls)
                value_matrix = self.__build_common_co_matrix__(cluster_lst, co_matrix, old_matrix=value_matrix)
                
            #end loop
            return cluster_lst
                
        finally:
            global score_dct
            score_dct.clear()
            
    def __split_bin__(self, cluster: Cluster) -> List[Cluster]:
        result = []
        for item in cluster:
            cluster2 = Cluster()
            cluster2.add(item)
            result.append(cluster2)
        return result
    
    def __refine_clusters__(self, source_lst: List[Cluster], value_matrix: SparseDictHashMatrix[Cluster, float],\
        skip_set: set) -> List[Cluster]:
        
        def get_score(cls: Cluster) -> float:
            global score_dct
            if cls not in score_dct: score_dct[cls] = self.bin_evaluator.score(cls) 
            return score_dct[cls]
        
        # match_lst = match_lst if match_lst is not None else lambda _: source_lst
        
        new_clusters = []
        for i in tqdm(range(len(source_lst))):
            c1 = source_lst[i]
            if c1 in skip_set or len(c1) == 0: continue
            score1 = get_score(c1)
            
            if score1 < 0.0:
                skip_set.add(c1)
                new_clusters += self.__split_bin__(c1)
                continue
            
            for j in range(i+1, len(source_lst)):
                c2 = source_lst[j]
                if c2 in skip_set or len(c2) == 0: continue
                sim = value_matrix.getEntry(c1, c2)
                if sim > self.min_threshold:
                    score2 = get_score(c2)
                    #Reminder that (c1 intersection c2) = Ø here,
                    #Scoring the item-chain will not result in polluted score.
                    combo = self.bin_evaluator.score_item_lst(itertools.chain(c1, c2))
                    if combo >= max(score1, score2):
                        skip_set.add(c1)
                        skip_set.add(c2)
                        mc = Cluster.merge(c1, c2)
                        new_clusters.append( mc )
                        break
        #end loop
        return new_clusters