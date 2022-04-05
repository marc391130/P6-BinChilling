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

    def read_file_fast(self, file_path: str, numpy_file: str) -> Dict[str, ContigData]:
        print("trying to load as numpy data...")
        numpy_data = self.try_load_numpy(numpy_file)
        if numpy_data is not None:
            print("found numpy file")
            return numpy_data
        else:
            print("could not find numpy file, loading from fasta...")
        data = self.read_file(file_path)
        print("saving fasta as numpy, so its faster in the future")
        self.save_numpy(data, numpy_file)
        return data  


    def read_file(self, file_path: str) -> Dict[str, ContigData]:
        result = {}
        abundance_length_dict = self.__get_abundance_length_dict__(const.ABUNDANCE_FILEPATH)

        def clean_line_name(line: str) -> str:
            return line.split('>')[1].replace('\n', '')
        

        composition_analyzer = CompositionAnalyzer()
        with open(file_path, 'r') as file:
            lines = file.readlines()

            for index in tqdm(range(len(lines))):
                line = lines[index]
                if not line.startswith('>'):
                    continue
                name = clean_line_name(line)
                composition = Composition()
                contig = ContigData(name, composition, 0, abundance_length_dict[name][0])
                temp_string = ""
                for i in range(index+1, len(lines)):
                    if lines[i].startswith('>'):
                        break
                    temp_string += lines[i]
                contig.contig_length = len(temp_string)
                composition_analyzer.analyze_composition(composition, temp_string)
                result[name] = contig
    
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

    def save_numpy(self, data: Dict[str, ContigData], outputfile: str) -> None:
        outputfile = outputfile if outputfile.endswith('.npy') else outputfile + ".npy"
        np.save(outputfile, data)
    

    def try_load_numpy(self, filename:str) -> Dict[str, ContigData] or None:
        try:
            return self.load_numpy(filename)
        except IOError as e:
            return None

    def load_numpy(self, filename:str) -> Dict[str, ContigData]:
        return np.load(filename, allow_pickle=True).item()