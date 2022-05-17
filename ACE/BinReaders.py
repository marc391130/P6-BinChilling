from os import listdir
from os.path import join
from EnsemblerTools import BinLogger
from ClusterDomain import Cluster, Partition, PartitionSet
from typing import Callable, Tuple, Dict, List
from Domain import ContigData, Composition
import Constants as const
from tqdm import tqdm
import itertools
import re
from multiprocessing import cpu_count, Pool
import numpy as np

class SCGReader:
    def __init__(self, filepaths: List[str] or str, DB_filepaths: List[str] = None, logger: BinLogger = None) -> None:
        if isinstance(filepaths, str): filepaths = [filepaths]
        self.filepaths = filepaths
        self.db_filepaths = DB_filepaths
        self.log = logger if logger is not None else print

    def read_contig_scg_superset(self) -> set:
        scgs = self.read_scg()
        return set(itertools.chain.from_iterable(scgs.values()))

    def read_MS_scgs(self) -> set:
        if self.db_filepaths is None or len(self.db_filepaths) == 0: return self.__read_SCG_db_set__(None)
        return set(itertools.chain.from_iterable((self.__read_SCG_db_set__(x) for x in self.db_filepaths) ))

    def read_scg(self) -> Dict[str, set]:
        result = {}
        if self.filepaths is None or len(self.filepaths) == 0:
            self.log("No SCG filepath supplied, skipping reading of SCGs, despite it being enabled")
            return dict()

        for file in self.filepaths:
            _temp = {}
            if file.endswith('.tsv'): _temp = self.__read_marker_gene_stat_file__(file)
            elif file.endswith('.scg'): _temp = self.__read_bacteria_scg_file__(file)

            for key, value in _temp.items():
                result[key] = result.get(key, set()).union(value)

        return result
        

    def __read_marker_gene_stat_file__(self, filepath) -> Dict[str, set]:
        self.log("Read Marker gene scg file!")
        def parse_SCG_from_line(contig_name: str, scg_line: str) -> set:
            scg_line = scg_line.replace('\n', '')
            
            result = set()
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
                    result.add(scg)
                
            return result
        
        
        result = {}
        with open(filepath, 'r') as f:
            for line in f.readlines():
                name = line.split('\t')[0]
                
                startindex, endindex = line.find('{'), line.rfind('}')
                SCG_str = line[startindex:endindex+1]
                
                result[name] = parse_SCG_from_line(name, SCG_str)
        return result

    def __read_bacteria_scg_file__(self, filepath) -> Dict[str, set]:
        temp: Dict[str, set] = {}
        self.log("Reading Bacteria file!")
        with open(filepath, 'r') as f:
            lines = f.readlines()

            for line in lines:
                split_line = line.split('\t')
                name, scg = split_line[0], split_line[1]
                contig_name = name.split('_')[0]
                scg_name = scg.replace('\n', '')

                _ = temp.get(contig_name, [])
                _.append(scg_name)
                temp[contig_name] = _

        result_dct = {item: set(value) for item, value in temp.items()}
        return result_dct

    def __read_SCG_db_set__(self, ms_file: str) -> set:
        print(self.db_filepaths, "Here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<")
        if ms_file is None or len(ms_file) == 0: return self.read_contig_scg_superset()
        string = ''
        result = set()
        with open(ms_file, 'r') as f:
            string = ''.join(f.readlines())
            result.update(re.findall("(?<=')([a-zA-Z0-9.,]+)(?=')", string))
        return result

class ContigFilter:
    def __init__(self, min_len: int) -> None:
        self.min_len = min_len
    def predicate(self, contig: ContigData) -> bool:
        return contig.contig_length > self.min_len


class ContigReader:
    def __init__(self, 
                 fasta_file: str, 
                 scg_reader: SCGReader,
                 depth_file: str = None, 
                 enable_analyse_contig_comp: bool = False,
                 numpy_file: str = None,
                 max_threads: int or None = None,
                 logger: BinLogger = None):
        self.fasta_file = fasta_file
        self.depth_file = depth_file
        self.numpy_file = numpy_file
        self.enable_analyse_contig_comp = enable_analyse_contig_comp
        self.max_threads = max_threads if max_threads is not None else cpu_count()
        self.SCG_reader = scg_reader
        self.log = logger if logger is not None else print

    def read_file_fast(self, numpy_file: str or None = None, load_SCGs:bool = False) -> Dict[str, ContigData]:
        numpy_path = numpy_file if numpy_file is not None else self.numpy_file
        
        self.log("Trying to load cache as numpy data...")
        if numpy_path is not None:        
            numpy_data = self.try_load_numpy(numpy_path)
            if numpy_data is not None:
                self.log("found cache file")
                return numpy_data
            else:
                self.log("could not find cache file, loading from fasta...")
        else:
            self.log('no cache file supplied, loading fasta. Contig objects will not be cached...')
        data = self.read_file(self.fasta_file, load_SCGs)
        if numpy_path is not None:
            self.log("saving contigs as numpy, so its faster in the future")
            self.save_numpy(data, numpy_path)
        else:
            self.log('skipping caching result')
        return data


    def read_file(self, file_path: str, load_SCGs:bool = False) -> Dict[str, ContigData]:
        abundance_length_dict = self.__get_abundance_length_dict__(self.depth_file) if self.depth_file is not None else None
        temp_result: Dict[str, List[str]] = {}
        temp_result_index_map : Dict[str, int] = {}
        
        def get_abundance(name: str) -> List[float]:
            nonlocal abundance_length_dict
            if abundance_length_dict is None: return 0.0
            return abundance_length_dict[name][0] if not self.depth_file.endswith('.npz')\
                else abundance_length_dict[str(temp_result_index_map[name])]
        
        def clean_line_name(line: str) -> str:
            return line.split('>')[1].replace('\n', '')
        
        self.log("Reading fasta file...")        
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
        self.log(f"Analysing contig compositions (using fast mode: {fast_mode_answer})")
        with Pool(min(self.max_threads, cpu_count())) as p:
            contig_lst: List[ContigData] = list(tqdm(p.imap_unordered(__partial_build_contig_multithread__, parameters), total=len(parameters)))
        result = {contig.name: contig for contig in contig_lst}
        
        if load_SCGs:
            self.log("loading SCGs...")
            SCG_dct = self.SCG_reader.read_scg()
            for contig_name, contig in tqdm(result.items()):
                if contig_name in SCG_dct:
                    contig.SCG_genes = set(SCG_dct[contig_name])
        else:
            self.log("skipping SCG load, as option disabled...")
        
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
            #print(data)
            result = list(data['arr_0'])
            data.close()
            return result

        result = {}
        if file_path.endswith('.npz'):
            npz_data = read_npz(file_path)
            
            for data_idx in range(len(npz_data)):
                result[str(data_idx)] = (np.median(npz_data[data_idx], axis=0), 1)
            
        else:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                
                for i in range(1, len(lines)):
                    data = lines[i].split('\t')
                    name = data[0]
                    length = int(float(data[1]))
                    abundance = list([float(data[3]), float(data[4])])
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

def __partial_build_contig_multithread__(tuple: Tuple[str, List[str], List[float], bool]) -> ContigData:
    name, dna_seq, abundance, enable_comp_analyse = tuple
    
    dna_string = ''.join(dna_seq).replace('\n', '')
    composition = Composition()
    if enable_comp_analyse:
        analyze_composition(dna_string, composition)
    
    contig = ContigData(name, composition, contig_length=len(dna_string), abundance=abundance)
    return contig


class ClusterReader:
    def __init__(self, file_path: str, contig_reader: ContigReader, numpy_file: str = None) -> None:
        self.contigReader = contig_reader
        self.numpy_file = numpy_file
        self.clusters = self.__read_clusters__(file_path)

    def __read_clusters__(self, file_path: str) -> List[Cluster]:
        result_clusters = []
        cluster_data_map = {}
        with open(file_path, 'r') as f:
            for line in tqdm(f.readlines()):
                split_line = line.split('\t')
                cluster_idx, edge_name = split_line[0], split_line[1].replace('\n', '')
                cluster_data_map[cluster_idx] = cluster_data_map[cluster_idx] + [edge_name] if cluster_idx in cluster_data_map else [edge_name]
        
        contig_scg_dct = self.contigReader.read_file_fast(self.numpy_file, True)
        for cluster_idx, edge_lst in cluster_data_map.items():
            cluster = Cluster()
            for edge in edge_lst:
                r = contig_scg_dct.get(edge, None)
                if r is not None: cluster.add(r)
            result_clusters.append(cluster)
        return result_clusters

class PartitionSetReader:
    def __init__(self, cluster_folder_path: str, contig_reader: ContigReader, file_predicate: Callable[[str], bool] = lambda x: True,\
            contig_filter: ContigFilter = None, logger: BinLogger = None) -> None:
        self.folder_path = cluster_folder_path
        self.file_predicate = file_predicate
        self.contig_reader = contig_reader
        self.contig_filter = contig_filter if contig_filter is not None else ContigFilter(0)
        self.logger = logger if logger is not None else print
        

    def read_partisionSet(self) -> PartitionSet:
        return self.read_file()

    def read_file(self, show_warnings: bool = True) -> PartitionSet:
        contig_dct = {name: contig for name, contig in self.contig_reader.read_file_fast(None, True).items() if self.contig_filter.predicate(contig)}
        files = [f for f in listdir(self.folder_path) if self.file_predicate(f)]
        
        partition_set = PartitionSet[ContigData](list(contig_dct.values()))
        
        warning_clusters = []
        print(f"reading {len(files)} partitions files...")
        for file in files:
            partition = partition_set.create_partition() #dont have to add it to partition after this
            self.__parse_partition__(file, partition, contig_dct)
            
            if len(partition) <= 1:
                warning_clusters.append(f"partition in file { file } only contains 1 cluster")
        
        if show_warnings:
            [print("> " +w) for w in warning_clusters]
        return partition_set

    def __parse_partition__(self, filename: str, partition: Partition, contigs: Dict[str, ContigData]) -> None:
        with open(join(self.folder_path, filename), 'r') as f:
            for line in tqdm(f.readlines()):
                cluster_name, edge_name = PartitionSetReader.__parse_cluster_line__(line) 
                if edge_name in contigs:
                    partition.add(cluster_name, contigs[edge_name])


    @staticmethod
    def __read_single_partition__(filename: str) -> Partition:
        partition = Partition()

        with open(filename, 'r') as f:
            for line in f.readlines():
                cluster_name, edge_name = PartitionSetReader.__parse_cluster_line__(line)

                partition.add(cluster_name, edge_name)
        return partition

                
    #returns tuple of cluster_name, edge_name
    @staticmethod
    def __parse_cluster_line__(line: str) -> Tuple[str, str]:
        split_line = line.split('\t')
        return (split_line[0], split_line[1].replace('\n', ''))
    
