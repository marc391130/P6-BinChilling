from os import listdir, getcwd
from os.path import join, exists
import string
from typing import Callable, List, Dict
import Constants as CONSTANT
import sys
from tqdm import tqdm
from ContigReader import ContigReader



class bin_writer:
    def __init__(self, fasta_file: str, cluster_file: str, output_path: str, file_predicate: Callable[[str], bool] = None, append_filename: str = None) -> None:
        self.fasta_file = fasta_file
        self.cluster_file = cluster_file
        self.file_predicate = file_predicate if file_predicate is not None else lambda x: True
        self.append_filename = append_filename if append_filename is not None else "bin_"
        self.output_path = output_path
        
    def __get_folder_content(self, folder_path: str) -> List[str]:
        return [f for f in listdir(self.input_folder_path) if self.file_predicate(f)]
    
    # key is contig name, value is contig string 
    def read_fasta(self) -> Dict[str, str]:
        current = ''
        result = {}
        
        def contig_cleaner(contingname: str) -> str:
            return contingname.replace('>', '').replace('\n', '')
        
        print("reading fasta...")
        
        with open(self.fasta_file, 'r') as f:
            for line in tqdm(f.readlines()):
                clean = contig_cleaner(line)
                if line.startswith('>') and clean not in result:
                    current = clean
                    result[current] = ''
                elif line.startswith('>') and clean in result:
                    raise Exception("this shouldnt be possible")
                else:
                    result[current] += line.replace('\n', '')
        return result
    
    #cluster name to all the clusters
    def read_clusters(self) -> Dict[str, List[str]]:
        
        result = {}
        print("reading clusters...")
        
        with open(self.cluster_file, 'r') as f:
            for line in tqdm(f.readlines()):
                splits = line.replace('\n', '').split('\t')
                
                if splits[0] not in result:
                    result[splits[0]] = []
                result[splits[0]] += [splits[1]]
        return result
    
    #key of values is binned contigname and value is dna string
    def write_fasta(self, outputpath: str, values: Dict[str, str]):
        n = 80 #how many chars per line
        print('writing new fasta file...')
        outputpath = outputpath if outputpath.endswith(".fasta") else outputpath + ".fasta"
        
        with open(outputpath, 'x') as f:
            for binname, dna_string in values.items():
                f.write('>' + binname + '\n')
                for dna_line in [dna_string[i:i+n] for i in range(0, len(dna_string), n)]:
                    f.write(dna_line + '\n')
    
    def write_bins(self, values: Dict[str, Dict[str, str]]):
        print("writing bins...")
        for clustername, cluster in tqdm(values.items()):
            name = self.append_filename + clustername + ".fasta"
            self.write_fasta(join(self.output_path, name), cluster)
        
    def work(self, showClusters: bool = True):
        
        clusters = self.read_clusters()
        contigs = self.read_fasta()
        
        #split contigs into clusters
        print("split contigs into clusters...")
        result = {binname: {name: contigs[name] for name in contignames} \
            for binname, contignames in tqdm(clusters.items())} #cluster name to dict of contigname to contingstring 
        
        self.write_bins(result)
        
        if showClusters:
            for k, v in result.items():
                print(f">{k}") 
            for k2, v2 in v.items():
                print(f"{k}_{k2} len: {len(v2)}")               
                
if __name__ == '__main__':
    print(sys.argv)
    writer = bin_writer(sys.argv[1], sys.argv[2], sys.argv[3], \
        lambda x: x.startswith('cluster'),
        sys.argv[4] if len(sys.argv) >= 5 else None)
    writer.work()