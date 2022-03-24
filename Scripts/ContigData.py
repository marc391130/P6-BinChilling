from unicodedata import name
from Composition import Composition
from typing import List

class ContigData: 
    def __init__(self, name: str = "", composition: Composition = None, contig_length: int = 0, abundance: float = 0):
        self.composition = composition
        self.name = name
        self.contig_length = contig_length
        self.abundance = abundance
    
    def as_composition_list(self, addatiive_value = 0) -> List[float]:
        result = []
        i = 0
        for key, value in self.composition.AsNormalized().items():
            result.append(value + value * addatiive_value * i)
            i += 1
        return result

    def pretty_print(self) -> None:
        print(f"{self.name} {self.contig_length} {self.abundance} {self.composition}")
    
    def __hash__(self) -> int:
        return self.name.__hash__()