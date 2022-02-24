from Composition import Composition

class ContigData: 
    def __init__(self, name: str = "", composition: Composition = None, contig_length: int = 0, abundance: float = 0):
        self.composition = composition
        self.name = name
        self.contig_length = contig_length
        self.abundance = abundance
    
    def pretty_print(self) -> None:
        print(f"{self.name} {self.contig_length} {self.abundance} {self.composition}")