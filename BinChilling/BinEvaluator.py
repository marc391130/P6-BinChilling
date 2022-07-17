import itertools
import sys
from math import log, sqrt
from typing import Iterable, Iterator, List, Dict, Tuple, TypeVar, Generic, Set
from tqdm import tqdm

from ClusterDomain import Cluster
from Domain import ContigData, bin_size
from SparseMatrix_implementations import SparseDictHashMatrix, SortKeysByHash
from EnsemblerTools import BinLogger
from Cluster_matrices import CoAssosiationMatrix
from scipy.sparse import csr_matrix, lil_matrix
from multiprocessing import Pool, cpu_count



class BinEvaluator:
    def __init__(self, all_SCGs: set  ) -> None:
        self.all_SCGs = all_SCGs
        
    def score(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> float:
        comp, conn, mp = self.evaluate(cluster, skip_item, include_item)
        return self.calc_score(comp, conn, mp)
    
    
        
    # def score_len(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> float:
    #     size = bin_size((x for x in self.__chain_cluster__(cluster, include_item=include_item) if x is not skip_item ))
    #     if len(self.genome_size_range) == 0: return 0.0
    #     return max( (self.__calc_len_score__( base_size - size ) for base_size in self.genome_size_range) )
        
    def __calc_len_score__(self, size: int) -> float:
        return 1 / (log(0.001 * abs(size) + 2, 2) )
    
    def score_SCG(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> float:
        comp, conn, mp = self.evaluate(cluster, skip_item, include_item)
        return self.calc_score(comp, conn, mp)
    
    def evaluate_lst(self, clusters: List[Cluster]) -> Dict[Cluster, Tuple[float, float, float]]:
        return {cluster: self.evaluate(cluster) for cluster in clusters } 
    
    def evaluate(self, cluster: Cluster, skip_item: ContigData = None,  include_item: ContigData = None) -> Tuple[float, float, float]:
        scg_count = self.__calculate_number_of_SCGs__(cluster, skip_item, include_item)
        return self.__evaluate_scg_count__(scg_count)
    
    
    def score_lst(self, cluster_lst: List[Cluster]) -> Dict[Cluster, float]:
        return {cluster: self.score(cluster) for cluster in cluster_lst}
    
    def score_item_lst(self, item_lst: Iterable[ContigData]) -> float:
        scg_count = self.__calculate_number_of_SCGs__(item_lst)
        return self.__score_scg_count__(scg_count)
    
    def calc_score(self, completeness: float, contamination: float, megabin_pen: float) -> float:
        return completeness - (contamination)**2 - sqrt(megabin_pen)
        # return completeness - 0.5*contamination - 0.5*megabin_pen
    
    def score_items(self, cluster: Cluster, extra_item: ContigData = None) -> Dict[ContigData, float]:
        if len(cluster) == 0 and extra_item is None: return {}
        if len(cluster) == 1 and (extra_item is None or extra_item in cluster) or\
            (len(cluster) == 0 and extra_item is not None): 
            return { item: self.score(cluster, include_item=extra_item) for item in self.__chain_cluster__(cluster, extra_item) }
        SCGs = self.__calculate_number_of_SCGs__(cluster, None, extra_item)
        total_value = self.__score_scg_count__(SCGs)

        result = {item: total_value -  self.__score_scg_count__(self.__remove_item_from_SCG_count__(item, SCGs))\
            for item in self.__chain_cluster__(cluster, extra_item)}
        return result
    
    def calculate_item_score(self, cluster: Cluster, extra_item: ContigData or None = None) -> Dict[ContigData, float]:
        if len(cluster) == 0 and extra_item is None: return {}
        if len(cluster) == 1 and (extra_item is None or extra_item in cluster) or\
            (len(cluster) == 0 and extra_item is not None): 
            return { item: self.score(cluster, include_item=extra_item) for item in self.__chain_cluster__(cluster, extra_item) }
        SCGs = self.__calculate_number_of_SCGs__(cluster, None, extra_item)
        total_value = self.__score_scg_count__(SCGs)

        result = {item: total_value-  self.__score_scg_count__(self.__remove_item_from_SCG_count__(item, SCGs))\
            for item in self.__chain_cluster__(cluster, extra_item)}
        return result
    
    
    
    def __remove_item_from_SCG_count__(self, item: ContigData, scgs: Dict[str, int]) -> Dict[int, str]:
        return { scg: ((count - 1) if scg in item.SCG_genes else count )\
            for scg, count in scgs.items() if not (count <= 1 and scg in item.SCG_genes)  }
    
    #?______________ PRIVATE METHODS ______________
    def __score_scg_count__(self, scg_count: Dict[str, int]) -> float:
        comp, conn, mp = self.__evaluate_scg_count__(scg_count)
        return self.calc_score(comp, conn, mp)
    
    def __evaluate_scg_count__(self, scg_count: Dict[str, int]) -> Tuple[float,float,float]:
        return \
            self.__calculate_completeness__(scg_count), \
            self.__calculate_contamination__(scg_count),\
            self.__calculate_megabin_penalty__(scg_count)
    
    def __calculate_completeness__(self, SCG_count: Dict[str, int]) -> float:
        
        counter = len(SCG_count)
        divisor = len(self.all_SCGs)
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_contamination__(self, SCG_count: Dict[str, int]) -> float:
        dSCG = [scg for scg, count in SCG_count.items() if count > 1]

        counter = len(dSCG)
        divisor = len(SCG_count) #This is the same as uniques, but does not require an additional n^2 call
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_purity__(self, SCG_count: Dict[str, int]) -> float:
        dSCG = [scg for scg, count in SCG_count.items() if count > 1]
        
        counter = len(SCG_count) #Same as uniques, without having to compute
        divisor = len(SCG_count) + len(dSCG)
        return (counter / divisor)*100 if divisor != 0 else 0

    def __calculate_megabin_penalty__(self, SCG_count: Dict[str, int]) -> float:
        nr_SCGs = sum(SCG_count.values()) - len(SCG_count)
        return (nr_SCGs / len(self.all_SCGs))*100 if len(self.all_SCGs) != 0 else 0
    
    def __calculate_number_of_SCGs__(self, cluster: Iterable[ContigData], skip_item: ContigData = None, include_item: ContigData = None) -> Dict[str, int]:
        SCGs = {}
        for item in self.__chain_cluster__(cluster, include_item):
            if skip_item is not None and item is skip_item: continue
            for scg in item.SCG_genes: 
                SCGs[scg] = SCGs.get(scg, 0) + 1
        return SCGs
    
    def __sum_size__(self, cluster: Cluster, skip_item: ContigData or None = None, include_item: ContigData or None = None):
        return sum((x.contig_length for x in self.__chain_cluster__(cluster, include_item) if x is not skip_item))
    
    def __chain_cluster__(self, cluster: Iterable[ContigData], include_item: ContigData or None = None) -> Iterator[ContigData]:
        return cluster.__iter__() if include_item is None or include_item in cluster\
            else itertools.chain(cluster.__iter__(), [include_item])
    


score_dct: Dict[Cluster, float] = {}
shared_co_dct: Dict[int, float] = {}

class BinRefiner:
    def __init__(self, bin_evaluator: BinEvaluator, min_threshold: float,
                 chunksize: int, logger: BinLogger) -> None:
        self.bin_evaluator = bin_evaluator
        self.log = logger
        self.min_threshold = min_threshold
        self.chunksize = chunksize
        self.score_dct: Dict[Cluster, float] = {}
    
    def __build_common_co_matrix__(self, cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix,\
        old_matrix: SparseDictHashMatrix[Cluster, Cluster, float] = None) -> SparseDictHashMatrix[Cluster, Cluster, float]:
        avg_co_matrix = SparseDictHashMatrix[Cluster, Cluster, float](SortKeysByHash, default_value=0.0)
        
        for i in tqdm(range(len(cluster_lst)), 'Adjusting CCO'):
            c1 = cluster_lst[i]
            for j in range(i+1, len(cluster_lst)):
                c2 = cluster_lst[j]
                value =  co_matrix.bin_mean(c1, c2) if old_matrix is None or old_matrix.has_entry(c1, c2) is False\
                    else old_matrix[c1, c2]
                avg_co_matrix.set_entry(c1, c2, value if value > self.min_threshold else 0.0)
        return avg_co_matrix
        
    def __build_shared_memory_co_matrix__(self, dictionary: Dict[int, float], co_matrix: CoAssosiationMatrix) \
        -> Tuple[Dict[int, Cluster], csr_matrix]:
        
        for item, index2_map  in tqdm(co_matrix.items()):
            # item, index2_map: Tuple[object, Dict[object, float]]
            for item2, value in index2_map.items():
                index = hash((hash(item), hash(item2)))
                if index in dictionary: raise Exception("non unique hash value optained from two items")
                dictionary[index] = value
            
    
    def __build_common_co_multiprocess__(self, cluster_lst: List[Cluster], new_clusters: List[Cluster], \
        matrix: SparseDictHashMatrix[Cluster, Cluster, float]) ->  SparseDictHashMatrix[Cluster, Cluster, float]:

        new_hash_lst = [[hash(item) for item in cluster] for cluster in new_clusters]
        
        if len(cluster_lst) > 0:
            old_hash_lst = [[hash(item) for item in cluster] for cluster in cluster_lst]
            parameters = ( ((i, j), new_hash_lst[i], old_hash_lst[j]) \
                for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P1') \
                    for j in range(len(cluster_lst)))
            
            #go through all new clusters and see if they are related to existing clusters
            with Pool(cpu_count()) as p:
                for id, value in p.imap_unordered(__partial_common_co_entry__, iterable=parameters, chunksize=self.chunksize):
                    id1, id2 = id
                    cluster1, cluster2 = new_clusters[id1], cluster_lst[id2]
                    matrix[cluster1, cluster2] = value
            del old_hash_lst, parameters
        #end if
        
        #new go through all the new clusters to see if they are related
        parameters2 = ( ((i, j), new_hash_lst[i], new_hash_lst[j]) \
            for i in tqdm(range(len(new_clusters)), 'adding clusters to CCO P2') \
                for j in range(i+1, len(new_clusters)))
        
        with Pool(cpu_count()) as p:
            for id, value in p.imap_unordered(__partial_common_co_entry__, iterable=parameters2, chunksize=self.chunksize):
                id1, id2 = id
                cluster1, cluster2 = new_clusters[id1], new_clusters[id2]
                matrix[cluster1, cluster2] = value      
        
        return matrix
                
    
    def refine_multiprocess(self, cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix) -> List[Cluster]:
        
        self.log(f'\n\nStarting bin refinement of {len(cluster_lst)} bins...')
        self.log(f'Building shared memory co_matrix...')
        global shared_co_dct
        self.__build_shared_memory_co_matrix__(shared_co_dct, co_matrix)
        
        self.log('Calculating cluster co-assosiation matrix using CO matrix...')
        value_matrix = self.__build_common_co_multiprocess__(\
            [], cluster_lst, SparseDictHashMatrix(SortKeysByHash, default_value= 0.0))
        
        try:
            while True:
                    remove_set = set()
                    new_clusters = self.__refine_clusters__(cluster_lst, value_matrix, remove_set)
                    if len(new_clusters) < 2:
                        break 
                    else:
                        self.log(f'Refined {len(remove_set)} clusters into {len(new_clusters)} new clusters.\n')
                        # self.log('Adjusting cluster co-assosiation matrix...')
                    
                    #remove skip set items
                    for cls in remove_set:
                        cluster_lst.remove(cls)
                    
                    value_matrix = self.__build_common_co_multiprocess__(cluster_lst, new_clusters, value_matrix)
                    
                    #add the new clusters to the cluster list
                    for cls in new_clusters:
                        cluster_lst.append(cls)
            #end while
            self.log('\n')
            return cluster_lst            
        finally:
            #reset shared memory
            shared_co_dct.clear()
        
        
    def Refine(self, cluster_lst: List[Cluster], co_matrix: CoAssosiationMatrix) -> List[Cluster]:
        
        self.log(f'\n\nStarting bin refinement of {len(cluster_lst)} bins...')
        self.log('Calculating cluster co-assosiation matrix using CO matrix...')
        value_matrix = self.__build_common_co_matrix__(cluster_lst, co_matrix)
        
        try:
            self.log('\nRefining bins...')
            source_clusters = cluster_lst    
            while True:
                skip_set = set()
                new_clusters = self.__refine_clusters__(cluster_lst, value_matrix, skip_set)
                if len(new_clusters) == 0:
                    break 
                else:
                    source_clusters = new_clusters
                    self.log(f'Refined {len(skip_set)} clusters into {len(new_clusters)} new clusters.\n')
                    # self.log('Adjusting cluster co-assosiation matrix...')
                
                for cls in skip_set:
                    cluster_lst.remove(cls)
                for cls in new_clusters:
                    cluster_lst.append(cls)
                value_matrix = self.__build_common_co_matrix__(cluster_lst, co_matrix, old_matrix=value_matrix)
                
            #end loop
            self.log('\n')
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
    
    def __refine_clusters__(self, source_lst: List[Cluster], value_matrix: SparseDictHashMatrix[Cluster, Cluster, float],\
        skip_set: set) -> List[Cluster]:
        
        def get_score(cls: Cluster) -> float:
            global score_dct
            if cls not in score_dct: score_dct[cls] = self.bin_evaluator.score(cls) 
            return score_dct[cls]
        
        # match_lst = match_lst if match_lst is not None else lambda _: source_lst
        
        new_clusters = []
        for i in tqdm(range(len(source_lst)), desc='Refining bins'):
            c1 = source_lst[i]
            if c1 in skip_set or len(c1) == 0: continue
            score1 = get_score(c1)
            
            for j in range(i+1, len(source_lst)):
                c2 = source_lst[j]
                if c2 in skip_set or len(c2) == 0: continue
                sim = value_matrix.getEntry(c1, c2)
                if sim <= self.min_threshold: continue
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
            
            if c1 in skip_set: continue
            if score1 < 0.0:
                skip_set.add(c1)
                new_clusters += self.__split_bin__(c1)
                continue
        #end loop
        return new_clusters
    
    
def __partial_common_co_entry__(tuple: Tuple[object, List[int], List[int]]) -> float:
    id, cluster1, cluster2 = tuple
    global shared_co_dct
    value = 0.0
    
    for item1 in cluster1:
        for item2 in cluster2:
            index = hash((item1, item2))
            if index in shared_co_dct:
                value += shared_co_dct[index]
    
    value = value / (len(cluster1) + len(cluster2))
    return (id, value)