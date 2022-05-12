from typing import Dict, List, Tuple
import itertools
import re

class SCGReader:
    def __init__(self, filepaths: List[str], DB_filepaths: List[str] = None) -> None:
        if isinstance(filepaths, str): filepaths = [filepaths]
        self.filepaths = filepaths
        self.db_filepaths = DB_filepaths
        self.db_scgs = None if DB_filepaths is None else set()

    def read_contig_scg_superset(self) -> set:
        scgs = self.read_scg()
        return set(itertools.chain(scgs.values()))

    def get_db_set_data(self):
        if self.db_scgs is None:
            raise Exception("You tried to get SCG DB data, however you have not supplied files or not read them yet using read_scg()!")
        return self.db_scgs

    def read_scg(self) -> Dict[str, set]:
        result = {}
        if self.filepaths is None or len(self.filepaths) == 0:
            print("No SCG filepath supplied, skipping reading of SCGs, despite it being enabled")
            return dict()

        for file in self.filepaths:
            _temp = {}
            if file.endswith('.tsv'): _temp = self.__read_marker_gene_stat_file__(file)
            elif file.endswith('.scg'): _temp = self.__read_bacteria_scg_file__(file)

            for key, value in _temp.items():
                result[key] = result.get(key, set()).union(value)

        if self.db_filepaths is None or len(self.db_filepaths) == 0:
            print("No DB SCG filepath supplied, skipping reading of DB SCGs")
        else:
            for file in self.db_filepaths:
                _temp = {}
                if file.endswith('.ms'): self.db_scgs = self.db_scgs.union(self.__read_SCG_db_set__(file))
                else: print(f"file: {file} is not a .ms file!")

        return result
        

    def __read_marker_gene_stat_file__(self, filepath) -> Dict[str, set]:
        print("Read Marker gene scg file!")
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
        print("Reading Bacteria file!")
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

    def __read_SCG_db_set__(self, file) -> set:
        string = ''
        result = set()
        for file in self.db_filepaths:
            with open(file, 'r') as f:
                string = ''.join(f.readlines())
            for item in set(re.findall("(?<=')([a-zA-Z0-9.,]+)(?=')", string)):
                result.add(item)
            string = ''

        return result

if __name__ == '__main__':
    scg_reader = SCGReader(['../Dataset/marker_gene_stats.tsv', '../Dataset/_proteins.faa.bacteria.scg'], ['../Dataset/Bacteria.ms']) # _proteins.faa.bacteria.scg
    data = scg_reader.read_scg()
    db_data = scg_reader.get_db_set_data()
    for item, value in data.items():
        print(item,value)

    # print(db_data)

    # for item in db_data:
    #     print(item)
