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
        with open(file_path, 'r') as file:
            lines = file.readlines()

            composition_analyzer = CompositionAnalyzer()
            temp_contig = ContigData()
            temp_composition = Composition()

            i = 0
            temp_string = ""
            for i in tqdm(range(len(lines))):
                if lines[i].startswith('>'):
                    if i != 0:
                        composition_analyzer.analyze_composition(temp_composition, temp_string)
                        temp_contig.composition = temp_composition
                        result[temp_contig.name] = temp_contig
                        temp_string = ""

                    name = lines[i].split('>')[1].replace("\n","")
                    temp_contig = ContigData(name, abundance= abundance_length_dict[name][0], contig_length=abundance_length_dict[name][1])
                    temp_composition = Composition()
                else:
                    temp_string = temp_string + lines[i]
        return result

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

    