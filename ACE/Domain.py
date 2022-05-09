import itertools
from math import floor
from typing import List, Dict, Tuple, Iterable
import Constants as const

class Composition(Dict[str, float]):
    def __init__(self):
        self.__recursive_setup__("", const.COMPOSITION_CONSTANT - 1)
        return
    
    def __recursive_setup__(self, parent_str : str, depth : int) -> str:
        for label in const.ACID_LABELS:
            key = parent_str + label
            
            if depth == 0:
                self[self.__sort_ACID__(key)] = 0
            else:
                self.__recursive_setup__(key, depth - 1)

    def __setItem__(self, key: str, value: float):
        self.__assert_ACID__(key)
        super().__setitem__(key, value)

    def __assert_ACID__(self, item: str) -> None:
        if len(item) != const.COMPOSITION_CONSTANT or const.ACID_SET.issuperset(item) is False:
            raise Exception(f"ACID string '{item}' doesnt match the ccomposition constant {const.COMPOSITION_CONSTANT}")

    def __sort_ACID__(self, key: str) -> str:
        reverse = key[::-1]
        return key if key < reverse else reverse

    def AddOccurence(self, key:str):
        k = self.__sort_ACID__(key)
        self.__assert_ACID__(k)
        self[k] += 1

    def AsNormalized(self) -> Dict[str, float]:
        total = sum([self[k] for k in self]) 
        dic = {} 
        for key, value in self.items():
            dic[key] = value / total
        return dic

    def AsNormalizedFeatureList(self) -> List[float]:
        total = sum([self[k] for k in self]) 
        return [x / total for x in self.values()]

       
    #prints only the values over the threshold
    def PrettyPrintValues(self, threshold: int = 0):
        print(f"printing compositions over threshold { threshold }")
        self.__pretty_print__(self, threshold)
    
    def PrettyPrintNormalized(self, threshold: float = 0):
        print(f"printing normalized compositions over threshold { threshold }")
        self.__pretty_print__(self.AsNormalized(), threshold)

    def __pretty_print__(self, composition: Dict[str, float], threshold:float):
        for key, value in composition.items():
            if threshold < value:
                print(f"composition of {key}: { value }")

class ContigData: 
    def __init__(self, name: str = "", composition: Composition = None, contig_length: int = 0, avg_abundance: float = 0):
        self.composition = composition
        self.name = name
        self.contig_length = contig_length
        self.avg_abundance = avg_abundance
        self.SCG_genes = set()
    
    def as_composition_list(self, addatiive_value = 0) -> List[float]:
        result = []
        i = 0
        for key, value in self.composition.AsNormalized().items():
            result.append(value + value * addatiive_value * i)
            i += 1
        return result

    def __has_analysis__(self) -> bool:
        return any( (x != 0.0 for x in self.composition.values()) )    

    def pretty_print(self) -> None:
        print(f"{self.name} {self.contig_length} {self.avg_abundance} {self.composition}")
    
    def __hash__(self) -> int:
        return self.name.__hash__()

def bin_size(contig_lst: Iterable[ContigData]) -> int:
    return sum( (x.contig_length for x in contig_lst) )
    
def compute_3_4_scg_count(contig_lst: List[ContigData]) -> int:
    return compute_3_4_avg(count_scgs(contig_lst)) 
    

def compute_3_4_avg(scg_count: Dict[str, int]) -> int:
    values = scg_count.values()
    avg = sum(values) / len(values)
    return floor(avg + ((max(values) - avg) / 2))


def count_scgs(contig_lst: List[ContigData]) -> Dict[str, int]:
    result = {}
    for scg in itertools.chain.from_iterable([x.SCG_genes for x in contig_lst]):
        result[scg] = result.get(scg, 0) + 1
    return result