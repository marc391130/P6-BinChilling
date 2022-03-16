from typing import Dict, List, Tuple
from Composition import Composition
from ContigData import ContigData
from tqdm import tqdm
from CompositionAnalyzer import CompositionAnalyzer
import Constants as const
import numpy as np

class ContigReader:
    def __init__(self):
        pass

    def read_file(self, file_path: str) -> Dict[str, ContigData]:
        result = {}
        abundance_length_dict = self.__get_abundance_length_dict__(const.ABUNDANCE_FILEPATH)

        def __reset_contig_composition__(name: str, abundance: float) -> Tuple[ContigData, Composition]:
            return (ContigData(name, abundance=abundance), Composition())

        def __handle_contig_string__(i: int, temp_contig: ContigData, temp_composition: Composition, temp_string: str, abundance_length_dict: Dict[str, Tuple[float, int]]) -> str:
            if i != 0:
                self.__assert_contig_length_equal__(abundance_length_dict[temp_contig.name][1], temp_string, temp_contig.name)
                composition_analyzer.analyze_composition(temp_composition, temp_string)
                temp_contig.composition = temp_composition
                temp_contig.contig_length = len(temp_string)
                result[temp_contig.name] = temp_contig
                return ""
            return temp_string

        with open(file_path, 'r') as file:
            lines = file.readlines()

            composition_analyzer = CompositionAnalyzer()
            temp_contig = ContigData()
            temp_composition = Composition()

            temp_string = ""
            for i in tqdm(range(len(lines))):
                if lines[i].startswith('>'):
                    name = lines[i].split('>')[1].strip('\n')
                    temp_string = temp_string.replace("\n", "")
                    
                    temp_string = __handle_contig_string__(i, temp_contig, temp_composition, temp_string, abundance_length_dict)
                    temp_contig, temp_composition = __reset_contig_composition__(name, abundance_length_dict[name][0])
                else:
                    temp_string = temp_string + lines[i]
        return result

    def read_contig_names(self, file_path: str) -> List[str]:
        with open(file_path, 'r') as file:
            return [line.split('>')[1].strip('\n') for line in file.readlines() if line.startswith('>')]


    
    def __assert_contig_length_equal__(self, depth_len:int, contig:str, name:str) -> None:

        if depth_len < 1000000:
            if depth_len != len(contig):
                raise Exception(f" On edge: {name} Depth len ({depth_len}) is not equal to string length ({len(contig)})")
        else:
            if (abs(depth_len - len(contig)) / depth_len) > 0.01:
                raise Exception(f"Over 1% error")

    def __get_abundance_length_dict__(self, file_path: str) -> Dict[str, Tuple[float,int]]:
        result = {}
        with open(file_path, 'r') as file:
            lines = file.readlines()
            
            for i in range(1, len(lines)):
                data = lines[i].replace('\t', " ").split(' ')
                name = data[0]
                length = int(float(data[1]))
                abundance = float(data[2])
                result[name] = (abundance, length)
        return result

    def __get_abundance__(self, abundance_dict: Dict[str, float], name: str) -> float:
        return abundance_dict[name]

    def save_as_numpy(self, fastafile:str, outputfile:str) -> None:
        outputfile = outputfile if outputfile.endswith('.npy') else outputfile + ".npy"
        np.save(outputfile, self.read_file(fastafile))

    def load_numpy(self, filename:str) -> Dict[str, ContigData]:
        return np.load(filename, allow_pickle=True).item()