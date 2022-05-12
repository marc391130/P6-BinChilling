from typing import Callable, Dict, List, Tuple
from Domain import ContigData, Composition
from tqdm import tqdm
import Constants as const
import re
import sys
import numpy as np
from multiprocessing import Pool, cpu_count
from SCGReader import SCGReader

class ContigFilter:
    def __init__(self, min_len: int) -> None:
        self.min_len = min_len
    def predicate(self, contig: ContigData) -> bool:
        return contig.contig_length > self.min_len


class ContigReader:
    def __init__(self, 
                 fasta_file: str, 
                 depth_file: str = None, 
                 SCG_filepath: List[str] = None,
                 SCG_db_path: List[str] = None,
                 enable_analyse_contig_comp: bool = False,
                 numpy_file: str = None,
                 max_threads: int or None = None):
        self.all_scg_db_path = SCG_db_path
        self.fasta_file = fasta_file
        self.SCG_filepath = SCG_filepath
        self.depth_file = depth_file
        self.numpy_file = numpy_file
        self.enable_analyse_contig_comp = enable_analyse_contig_comp
        self.max_threads = max_threads if max_threads is not None else cpu_count()
        self.SCG_reader = SCGReader(SCG_filepath, SCG_db_path)

    def read_file_fast(self, numpy_file: str or None = None, load_SCGs:bool = False) -> Dict[str, ContigData]:
        numpy_path = numpy_file if numpy_file is not None else self.numpy_file
        
        print("Trying to load cache as numpy data...")
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

        parameters = [(contig_name, dna_parts, get_abundance(contig_name), self.enable_analyse_contig_comp) \
            for contig_name, dna_parts in temp_result.items()]
        contig_lst = []
        fast_mode_answer = 'no' if self.enable_analyse_contig_comp else 'yes'
        print(f"Analysing contig compositions (using fast mode: {fast_mode_answer})")
        with Pool(min(self.max_threads, cpu_count())) as p:
            contig_lst: List[ContigData] = list(tqdm(p.imap_unordered(__partial_build_contig_multithread__, parameters), total=len(parameters)))
        result = {contig.name: contig for contig in contig_lst}
        
        if load_SCGs:
            print("loading SCGs...")
            SCG_dct = self.SCG_reader.read_scg()
            print("I do this!")
            for contig_name, contig in tqdm(result.items()):
                if contig_name in SCG_dct:
                    contig.SCG_genes = set(SCG_dct[contig_name])
        else:
            print("skipping SCG load, as option disabled...")
        
        return result
        

    # def read_contig_SCGs(self) -> Dict[str, set]:
    #     if self.SCG_filepath is None or len(self.SCG_filepath) == 0:
    #         print("No SCG filepath supplied, skipping reading of SCGs, despite it being enabled")
    #         return dict()

    #     if len(self.SCG_filepath) > 1:
    #         result = {}
    #         for file in self.SCG_filepath:
    #             data = SCGReader(file).read_scg()
    #             for edge, scg_set in data.items():
    #                 result[edge] = set([scg for scg in result.get(edge, [])] + [scg for scg in scg_set])

    #         return result

    #     self.SCG_filepath = self.SCG_filepath[0]
        
        

    # def read_total_SCGs_set(self) -> set:
    #     string = ''
    #     result = set()
    #     for file in self.all_scg_db_path:
    #         with open(file, 'r') as f:
    #             string = ''.join(f.readlines())
    #         for item in set(re.findall("(?<=')([a-zA-Z0-9.,]+)(?=')", string)):
    #             result.add(item)
    #         string = ''

    #     return result
        
    
    
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


def analyze_composition(contig_string : str, composition_obj = Composition()) -> None:
    com_len = len(contig_string) - (const.COMPOSITION_CONSTANT - 1)
    #com_normalize_val = 1 / com_len
    if com_len < const.COMPOSITION_CONSTANT:
        raise Exception("composition string is smaller than the COMPOSITION_CONSTANT: " + str(const.COMPOSITION_CONSTANT))

    for x in range(com_len):
        key = contig_string[x:x + const.COMPOSITION_CONSTANT]
        reversed_key = key[::-1]
        if key in composition_obj:
            composition_obj.AddOccurence(key)
        elif reversed_key in composition_obj:
            composition_obj.AddOccurence(reversed_key)

def __partial_build_contig_multithread__(tuple: Tuple[str, List[str], float, bool]) -> ContigData:
    name, dna_seq, abundance, enable_comp_analyse = tuple
    
    dna_string = ''.join(dna_seq).replace('\n', '')
    composition = Composition()
    if enable_comp_analyse:
        analyze_composition(dna_string, composition)
    
    contig = ContigData(name, composition, contig_length=len(dna_string), avg_abundance=abundance)
    return contig

if __name__ == "__main__":
    reader = ContigReader(sys.argv[1], sys.argv[2], None)
    r = reader.read_file(sys.argv[1], True)
    with open('./output.tsv') as f:
        for name, contig in r.items():
            f.writelines(f'{name}\t{contig.contig_length}')
    print("DONE")
    print(len(r))