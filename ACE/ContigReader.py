from typing import Callable, Dict, List, Tuple
from Domain import ContigData, Composition
from tqdm import tqdm
import Constants as const
import re
import numpy as np
from multiprocessing import Pool, cpu_count

CONTIG_MIN_LEN = 1000
class ContigFilter:
    def __init__(self, min_len: int) -> None:
        self.min_len = min_len
    def predicate(self, contig: ContigData) -> bool:
        return contig.contig_length > self.min_len

class CompositionAnalyzer:
    def __init__(self) -> None:
        pass

    def analyze_composition(self, composition_dict : Composition, contig_string : str) -> None:
        com_len = len(contig_string) - (const.COMPOSITION_CONSTANT - 1)
        #com_normalize_val = 1 / com_len
        if com_len < const.COMPOSITION_CONSTANT:
            raise Exception("composition string is smaller than the COMPOSITION_CONSTANT: " + str(const.COMPOSITION_CONSTANT))

        for x in range(com_len):
            key = contig_string[x:x + const.COMPOSITION_CONSTANT]
            reversed_key = key[::-1]
            if key in composition_dict:
                composition_dict.AddOccurence(key)
            elif reversed_key in composition_dict:
                composition_dict.AddOccurence(reversed_key)

class ContigReader:
    def __init__(self, 
                 fasta_file: str, 
                 depth_file: str = None, 
                 SCG_filepath: str = None, 
                 numpy_file: str = None,
                 max_threads: int or None = None, 
                 contig_filter: ContigFilter = None ):
        self.fasta_file = fasta_file
        self.SCG_filepath = SCG_filepath
        self.depth_file = depth_file
        self.numpy_file = numpy_file
        self.max_threads = max_threads if max_threads is not None else cpu_count()
        self.contig_filter = contig_filter if contig_filter is not None else ContigFilter(CONTIG_MIN_LEN)

    def read_file_fast(self, numpy_file: str or None = None, load_SCGs:bool = False) -> Dict[str, ContigData]:
        numpy_path = numpy_file if numpy_file is not None else self.numpy_file
        
        print("trying to load cache as numpy data...")
        if numpy_path is not None:        
            numpy_data = self.try_load_numpy(numpy_path)
            if numpy_data is not None:
                print("found cache file")
                return numpy_data
            else:
                print("could not find cache file, loading from fasta...")
        else:
            print('no cache file supplied, loading fasta. Contig objects will not be cached...')
        data = self.read_file(self.fasta_file, load_SCGs)
        if numpy_path is not None:
            print("saving contigs as numpy, so its faster in the future")
            self.save_numpy(data, numpy_path)
        else:
            print('skipping caching result')
        return data


    def read_file(self, file_path: str, load_SCGs:bool = False) -> Dict[str, ContigData]:
        return self.read_file_multithreaded(file_path, load_SCGs)
        result: Dict[str, ContigData] = {}
        abundance_length_dict = self.__get_abundance_length_dict__(self.depth_file) if self.depth_file is not None else None

        def clean_line_name(line: str) -> str:
            return line.split('>')[1].replace('\n', '')
        
        composition_analyzer = CompositionAnalyzer()
        with open(file_path, 'r') as file:
            lines = file.readlines()
            current_contig = 0
            for index in tqdm(range(len(lines))):
                line = lines[index]
                if not line.startswith('>'):
                    continue
                name = clean_line_name(line) 
                composition = Composition()
                abundance = abundance_length_dict[name][0] if not self.depth_file.endswith('.npz') else abundance_length_dict[str(current_contig)][0]
                contig = ContigData(name, composition, 0, abundance)
                temp_string = ""
                current_contig += 1
                for i in range(index+1, len(lines)):
                    if lines[i].startswith('>'):
                        break
                    temp_string += lines[i]
                contig.contig_length = len(temp_string)
                composition_analyzer.analyze_composition(composition, temp_string)
                result[name] = contig
        if load_SCGs:
            print("loading SCGs...")
            SCG_dct = self.read_contig_SCGs()
            for contig_name, contig in tqdm(result.items()):
                if contig_name in SCG_dct:
                    contig.SCG_genes = set(SCG_dct[contig_name])
        else:
            print("skipping SCG load, as option disabled...")

        return result

    def read_file_multithreaded(self, file_path: str, load_SCGs:bool = False) -> Dict[str, ContigData]:
        abundance_length_dict = self.__get_abundance_length_dict__(self.depth_file) if self.depth_file is not None else None
        temp_result: Dict[str, List[str]] = {}
        temp_result_index_map : Dict[str, int] = {}
        
        def get_abundance(name: str) -> float:
            return abundance_length_dict[name][0] if not self.depth_file.endswith('.npz')\
                else abundance_length_dict[str(temp_result_index_map[name])][0]
        
        def clean_line_name(line: str) -> str:
            return line.split('>')[1].replace('\n', '')
        
        print("Reading fasta file...")        
        with open(file_path, 'r') as file:
            lines = file.readlines()
            current_name, current_index = '', 0
            for line in tqdm(lines):
                if line.startswith('>'):
                    current_name = clean_line_name(line)
                    temp_result[current_name] = []
                    temp_result_index_map[current_name] = current_index
                    current_index += 1
                else:
                    temp_result[current_name].append(line)

        parameters = [(key, value, get_abundance(key)) for key, value in temp_result.items()]
        contig_lst = []
        print("Analysing contig compositions...")
        with Pool(min(self.max_threads, cpu_count())) as p:
            contig_lst: List[ContigData] = [x for x in\
                tqdm(p.imap(__build_contig_multithread__, parameters), total=len(parameters)) if self.contig_filter.predicate(x)]
        result = {contig.name: contig for contig in contig_lst}
        
        if load_SCGs:
            print("loading SCGs...")
            SCG_dct = self.read_contig_SCGs()
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
                    raise Exception("loop count broke")
                if start == end or start+1 == end:
                    continue

                scg = scg_line[start+1:end]
                if scg.startswith(contig_name) is False:
                    result.append(scg)
                
            return result
        
        
        result = {}
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
            result = result.union(SCG_set)
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
        def read_npz(file_path):
            data = np.load(file_path)
            result = data['arr_0']
            data.close()
            return result

        result = {}
        if file_path.endswith('.npz'):
            npz_data = read_npz(file_path)
            
            for data_idx in range(len(npz_data)):
                result[str(data_idx)] = (np.median(npz_data[data_idx], axis=0), 0)
            
        else:
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
        wrapper = DataWrapper(self.fasta_file, data)
        outputfile = outputfile if outputfile.endswith('.npy') else outputfile + ".npy"
        np.save(outputfile, np.array([wrapper]))
    

    def try_load_numpy(self, filename:str) -> Dict[str, ContigData] or None:
        try:
            wrapper: DataWrapper = self.load_numpy(filename)[0]
            return wrapper.validate_and_get_data(self.fasta_file)
        except IOError as e:
            return None

    def load_numpy(self, filename:str) -> Dict[str, ContigData]:
        return np.load(filename, allow_pickle=True)
    
class DataWrapper():
    def __init__(self, identifier: str, data: np.ndarray):
        self.identifier = DataWrapper.cleanse_identifier(identifier)
        self.data = data

    def validate_and_get_data(self, identifier):
        cleansed_identifier = DataWrapper.cleanse_identifier(identifier)
        if cleansed_identifier != self.identifier:
            raise Exception(f"Cache identifier does not match (.fasta /.fna) file.\n \
                            Expected to load cache for {self.identifier}, but recieved {cleansed_identifier}.\n \
                            Please give another cache path or delete the old cache.")
        
        return self.data
    
    @staticmethod
    def cleanse_identifier(identifier :str) -> str:
        split_filename = identifier.split("/")
        split_filename = split_filename[len(split_filename) - 1].split("\\")
        return split_filename[len(split_filename) - 1]


def __build_contig_multithread__(tuple: Tuple[str, List[str], float]) -> ContigData:
    name, dna_seq, abundance = tuple
    
    dna_string = ''.join(dna_seq).replace('\n', '')
    
    composition = Composition()
    # analyser = CompositionAnalyzer()
    # analyser.analyze_composition(composition, dna_string)
    contig = ContigData(name, composition, contig_length=len(dna_string), avg_abundance=abundance)
    return contig

if __name__ == "__main__":
    reader = ContigReader('../Dataset/edges.fasta', '../Dataset/edges_depth.txt', '../Dataset/marker_gene_stats.tsv')
    r = reader.read_file_multithreaded('../Dataset/edges.fasta', True)
    print("DONE")
    print(len(r))