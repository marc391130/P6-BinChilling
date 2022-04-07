from typing import Callable, List
import sys
from tqdm import tqdm


COMPLETENESS_INDEX = 10
CONTAMINATION_INDEX = 11 
DEFAULT_PREDICATE = lambda x: True
BETA_MULTIPLIER = 5
QUALITY_THREDSHOLD = 0.5

class BinDto:
    def __init__(self, name: str, contamination: float, completeness: float) -> None:
        self.name = name
        self.contamination = contamination
        self.completeness = completeness
        
    def isNear(self) -> bool:
        return self.contamination <= 5 and self.completeness >= 90
    
    def isSubstantial(self) -> bool:
        return self.contamination > 5 and self.contamination <= 10 and\
            self.completeness < 90 and self.completeness >= 70
    
    def isModerate(self) -> bool:
        return self.contamination > 10 and self.contamination <= 15 and\
            self.completeness < 70 and self.completeness >= 50
    
    def isPartial(self) -> bool:
        return self.contamination > 15 and self.completeness < 50

    def checkQuality(self) -> float:
        quality = self.completeness / (self.contamination * BETA_MULTIPLIER + self.completeness) if self.completeness != 0 or self.contamination != 0 else 0
        return quality

    def isNotZeroCompleteness(self) -> bool:
        return self.completeness == 0
        
    def categoryToString(self) -> str:
        if self.isNear():
            return 'near'
        if self.isSubstantial():
            return 'substantial'
        if self.isModerate():
            return 'moderate'
        if self.isPartial():
            return 'partial'
        return 'Bad bin'
        
    def __str__(self) -> str:
        return f"{self.name}\t{self.contamination}\t{self.completeness}\t{self.categoryToString()}\t{self.checkQuality()}"
    
    @staticmethod
    def toStringFormat():
        return "name\tcontamination\tcompleteness\tcategory\tQuality\t"

class CheckMFilter:
    def __init__(self, filepath: str, outputPath: str, print_predicate: Callable[[BinDto], bool] = None) -> None:
        self.filepath = filepath
        self.outputPath = outputPath
        self.predicate = print_predicate if print_predicate is not None else DEFAULT_PREDICATE
    
    def run(self) -> List[BinDto]:
        result = []
        print("reading input...")
        with open(self.filepath, 'r') as f:
            for line in tqdm(f):
                if line.startswith('bin') is False:
                    continue
                splits = line.split(',')
                name = line.split('\t')[0]
                contamination = float(splits[CONTAMINATION_INDEX].split(":")[1])
                completeness = float(splits[COMPLETENESS_INDEX].split(":")[1])
                result.append(BinDto(name, contamination, completeness))
                
        return result
    
    def write_output(self, values: List[BinDto]) -> None:
        print("writing output...")
        with open(self.outputPath, 'w') as f:
            f.write(BinDto.toStringFormat() + "\n")
            for dto in tqdm(values):
                if self.predicate(dto):
                    f.write(str(dto) + "\n")
        
        print("Near:", len([dto for dto in values if dto.isNear()]))
        print("Substantial:", len([dto for dto in values if dto.isSubstantial()]))
        print("Moderate:", len([dto for dto in values if dto.isModerate()]))
        print("Partial:", len([dto for dto in values if dto.isPartial()]))
        print("Rest:", len([dto for dto in values if getPrintPredicate('rest')(dto)]))

            
    
    def prettyprint(self, line: str):
        splitlist = line.split(',')
        for split_index in range(len(splitlist)):
            split = splitlist[split_index]
            print("> " + split + "[" + str(split_index) + "]")
            
        print("________________|newline|_________________")



def getPrintPredicate(input: str = None) -> Callable[[BinDto], bool]:
    if input is None:
        return DEFAULT_PREDICATE
    elif input == "near":
        return lambda x: x.isNear()
    elif input == "substantial":
        return lambda x: x.isSubstantial()
    elif input == "moderate":
        return lambda x: x.isModerate()
    elif input == "partial":
        return lambda x: x.isPartial()
    elif input == "moderateOrBetter":
        return lambda x: x.isModerate() or x.isSubstantial() or x.isNear()
    elif input == "rest":
        return lambda x: x.isNear() is False and x.isModerate() is False and x.isSubstantial() is False and x.isPartial() is False and x.isNotZeroCompleteness() is False
    elif input == "all":
        return lambda x: x.isNear() or x.isModerate() or x.isSubstantial() or x.isPartial()
    elif input == "quality":
        return lambda x: x.checkQuality() >= QUALITY_THREDSHOLD
    else:
        raise Exception(f"print predicate '{input}' cannot be found."\
            + " try 'near', 'substantial', 'moderate', 'partial', 'moderateOrBetter' or leave blank")

if __name__ == "__main__":
    print(sys.argv)
    prediacte = getPrintPredicate(sys.argv[3]) if len(sys.argv) >= 4 else DEFAULT_PREDICATE
    writer = CheckMFilter(sys.argv[1], sys.argv[2], prediacte)
    writer.write_output(writer.run())