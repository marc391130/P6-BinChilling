import itertools
from typing import List, Dict, Tuple, TypeVar, Generic, Set
from tqdm import tqdm
from ClusterDomain import Cluster, Partition
from SparseMatrix_implementations import SparseDictHashMatrix, SparseTupleHashMatrix, SortKeysByHash, HashedMatrix
from EnsemblerTools import BinLogger, print_result, build_partition
from Cluster_matrices import CoAssosiationMatrix
import CoAssosiationFunctions as CoFunctions
from BinEvaluator import BinEvaluator
import os.path as path
from SharedDatastructures_implementations import HashTableIterator, SharedHashTable
import subprocess
import os

class RefinerBase:
    def __init__(self, bin_evaluator: BinEvaluator, output_filepath: str) -> None:
        self.output_filepath = output_filepath
    
    def refine(self, partition: Partition) -> None:
        raise NotImplementedError("You need to override this method")
        

class BinRefiner(RefinerBase):
    def __init__(self, bin_evaluator: BinEvaluator, min_threshold: float,
                 chunksize: int, logger: BinLogger) -> None:
        self.bin_evaluator = bin_evaluator
        self.log = logger
        self.min_threshold = min_threshold
        self.chunksize = chunksize
        self.score_dct: Dict[Cluster, float] = {}
    

    def refine(self, partition: Partition) -> Partition:
        cluster_lst = list(partition.values())
        self.log(f'\n\nStarting bin refinement of {len(cluster_lst)} bins...')
        
        self.log('Calculating cluster co-assosiation matrix using CO matrix...')
        co_cache = SparseTupleHashMatrix(SortKeysByHash, default_value=0)

        print(f"totalLen: {sum( (len(item.SCG_genes) for cluster in partition.values() for item in cluster ) )}")

        try:
            while True:
                remove_set = set()
                new_clusters = self.__refine_clusters_cacheless__(cluster_lst, co_cache, remove_set)
                if len(new_clusters) < 2:
                    break 
                else:
                    self.log(f'Refined {len(remove_set)} clusters into {len(new_clusters)} new clusters.\n')
                    # self.log('Adjusting cluster co-assosiation matrix...')
                
                #remove skip set items
                self.log('Releasing unkept memory')
                for cls in remove_set:
                    cluster_lst.remove(cls)
                    self.score_dct.pop(cls, 'default') #default such that an error is not thrown
                
                co_cache.pop_set(remove_set)
                
                #add the new clusters to the cluster list
                for cls in new_clusters:
                    cluster_lst.append(cls)
            #end while
            self.log('\n')
            return cluster_lst
        finally:
            pass
        
    def __split_bin__(self, cluster: Cluster) -> List[Cluster]:
        result = []
        for item in cluster:
            cluster2 = Cluster()
            cluster2.add(item)
            result.append(cluster2)
        return result
    
    
    def __refine_clusters_cacheless__(self, source_lst: List[Cluster],\
        co_cache: SparseTupleHashMatrix[Cluster, float],\
        skip_set: set) -> List[Cluster]:
        
        def get_score(cls: Cluster) -> float:
            if cls not in self.score_dct: 
                self.score_dct[cls] = self.bin_evaluator.score(cls) 
            return self.score_dct[cls]

        # match_lst = match_lst if match_lst is not None else lambda _: source_lst
        
        new_clusters = []
        for i in tqdm(range(len(source_lst)), desc='Refining bins'):
            c1 = source_lst[i]
            if c1 in skip_set or len(c1) == 0: continue
            score1 = get_score(c1)
            co_values = CoFunctions.Get_common_co_multiprocess(c1, 
                                                   source_lst[i+1:len(source_lst)], 
                                                   co_cache,
                                                   skip_set)
            
            for c2, sim in co_values.items():
                # if c2 in skip_set or len(c2) == 0: continue #This should be handled by the CoFunctions function...
                if sim <= self.min_threshold: continue
                score2 = get_score(c2)
                #Reminder that (c1 intersection c2) = Ø here,
                #Scoring the item-chain will not result in polluted score.
                combo = self.bin_evaluator.score_item_lst(itertools.chain(c1, c2))
                if combo > max(score1, score2):
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
            
            #save simularity
            co_cache.set_dict(c1, co_values)
            
        #end loop
        return new_clusters
    
    
EXECUTABLE = path.join(os.getcwd(), 'BinChillingTools.exe') 
class ExternalBinRefiner(RefinerBase):
    def __init__(self, partition_count: str,
                 output_filepath: str, 
                 partition_path: str, 
                 co_cache_path: str, 
                 logger: BinLogger) -> None:
        self._k = partition_count
        self._output_filepath = output_filepath
        self._co_cache_path = co_cache_path
        self._partition_path = partition_path
        self.log = logger
        self._connected = False
    
    def connect_module(self, throw_on_err: bool = True) -> None:
        out = subprocess.run([EXECUTABLE, '__test__'], check=True, text=True, capture_output=True)
        msg = out.stdout.split('\n')
        self.log(msg[0] if len(msg) > 0 else "Cannot connect to refiner module")
        success = '__REFINER_MODULE_CONNECTED__' in out.stdout
        
        self._connected = success
        if self._connected is False and throw_on_err:
            raise Exception('CANNOT CONENCT TO REFINER')
    
    def refine(self, partition: Partition) -> None:
        co_filename = CoFunctions.tmp_co_filename
        
        if self._connected is False:
            self.log("Refiner is not connected. Refinement step is skipped")
            return
        
        scg_filename = co_filename + '.scg'
        with open(scg_filename, 'w') as f:
            for item in (item for cluster in partition.values() for item in cluster):
                for scg in item.SCG_genes:
                    f.write( f"{item.name}\t{scg}\n")
            f.flush()
        subprocess.run([EXECUTABLE,
                        self._partition_path,
                        co_filename,
                        scg_filename,
                        str(float(self._k)),
                        self._output_filepath])
        try:
            os.remove(co_filename)
            os.remove(scg_filename)
        except:
            print("An error occured while trying to remove cache files.")
        
        
                    
            
        
        
        
        
        
        
        
        