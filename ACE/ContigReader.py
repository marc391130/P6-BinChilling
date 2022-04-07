from typing import Dict, List, Tuple
from Composition import Composition
from ContigData import ContigData
from tqdm import tqdm
from CompositionAnalyzer import CompositionAnalyzer
import re
import numpy as np

class ContigReader:
    def __init__(self, fasta_file: str, depth_file: str,  SCG_filepath: str = None, numpy_file: str = None):
        self.fasta_file = fasta_file
        self.SCG_filepath = SCG_filepath
        self.depth_file = depth_file
        self.numpy_file = numpy_file
        pass

    def read_file_fast(self, numpy_file: str = None, load_SCGs:bool = False) -> Dict[str, ContigData]:
        numpy_path = numpy_file if numpy_file is not None else self.numpy_file
        
        print("trying to load as numpy data...")
        if numpy_path is not None:        
            numpy_data = self.try_load_numpy(numpy_path)
            if numpy_data is not None:
                print("found npz file")
                return numpy_data
            else:
                print("could not find npz file, loading from fasta...")
        else:
            print('no npz path supplied, loading fasta. Contig objects will not be cached...')
        data = self.read_file(self.fasta_file, load_SCGs)
        if numpy_path is not None:
            print("saving fasta as numpy, so its faster in the future")
            self.save_numpy(data, numpy_path)
        else:
            print('skipping caching result')
        return data


    def read_file(self, file_path: str, load_SCGs:bool = False) -> Dict[str, ContigData]:
        result: Dict[str, ContigData] = {}
        abundance_length_dict = self.__get_abundance_length_dict__(self.depth_file)

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
        if load_SCGs:
            print("loading SCGs...")
            SCG_dct = self.read_SCGs()
            for contig_name, contig in tqdm(result.items()):
                if contig_name in SCG_dct:
                    contig.SCG_genes = set(SCG_dct[contig_name])
        else:
            print("skipping SCG load, as option disabled...")

        return result


    def read_contig_SCGs(self) -> Dict[str, List[str]]:
        if self.SCG_filepath is None:
            print("No SCG filepath supplied, skipping reading of SCGs, despite it being enabled")
            return dict()
        
        def parse_SCG_from_line(contig_name: str, scg_line: str) -> List[str]:
            scg_line = scg_line.replace('\n', '')
            
            result = []
            start_indecies = [int(i.start()) for i in re.finditer('\'', scg_line)]
            end_indecies = [int(i.start()) for i in re.finditer('\'', scg_line)]
            
            if len(start_indecies) != len(end_indecies):
                raise ValueError(f"{len(start_indecies)} != {len(end_indecies)}")
                        
            for i in range(0, len(start_indecies), 2):
                start, end = start_indecies[i], end_indecies[i+1]
                if start > end:
                    raise Exception("Kill me!")
                if start == end or start+1 == end:
                    continue

                scg = scg_line[start+1:end]
                if scg.startswith(contig_name) is False:
                    result.append(scg)
                
            return result
        
        
        result = {}
        print("Reading SCGs...")
        with open(self.SCG_filepath, 'r') as f:
            for line in f.readlines():
                name = line.split('\t')[0]
                
                startindex, endindex = line.find('{'), line.rfind('}')
                SCG_str = line[startindex:endindex+1]
                
                result[name] = parse_SCG_from_line(name, SCG_str)
        return result



    def read_total_SCGs_set(self) -> set:
        contig_SCGs = self.read_contig_SCGs()
        result = set()
        
        for SCG_lst in contig_SCGs.values():
            SCG_set = set(SCG_lst)
            result.union(SCG_set)
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
    
    
if __name__ == "__main__":
    reader = ContigReader('../Dataset/edges.fasta', '../Dataset/edges_depth.txt', '../Dataset/marker_gene_stats.tsv')
    r = reader.read_SCGs()
    print("DONE")
    for n, l in sorted(r.items(), key=lambda b: b[0]):
        print(n, '|', l)