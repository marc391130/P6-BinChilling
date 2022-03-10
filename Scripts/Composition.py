import Constants as const
from typing import Dict, List, TypeVar

class Composition(Dict[str, float]):
    def __init__(self):
        self.__recursive_setup__("", const.COMPOSITION_CONSTANT - 1)
        return
    
    def __recursive_setup__(self, parent_str : str, depth : int) -> str:
        for label in const.ACID_LABELS:
            key = parent_str + label
            reversed_key = key[::-1]
            
            if depth == 0:
                if reversed_key not in self:
                    self[key] = 0
            else:
                self.__recursive_setup__(key, depth - 1)
    

    def __setItem__(self, key: str, value: float):
        self.__assert_ACID__(key)
        super().__setitem__(key, value)

    def __assert_ACID__(self, item: str) -> None:
        cleanedItem = item
        for label in const.ACID_LABELS:
                cleanedItem = cleanedItem.replace(label, '')
        if len(item) != const.COMPOSITION_CONSTANT or len(cleanedItem) != 0:
            raise Exception(f"The item '{item}' contains other characters than ACID_LABELS: {const.ACID_LABELS}")

    def AddOccurence(self, key:str):
        self.__assert_ACID__(key)
        self[key] += 1

    def AsNormalized(self) -> Dict[str, float]:
        total = sum([self[k] for k in self]) 
        dic = {} 
        for key, value in self.items():
            dic[key] = value / total
        return dic

       
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
