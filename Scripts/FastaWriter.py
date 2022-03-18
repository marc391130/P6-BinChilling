from os import listdir, getcwd
from os.path import join, exists
import string
from typing import Callable, List, Dict
import Constants as CONSTANT
import sys
from tqdm import tqdm
from ContigReader import ContigReader



class fasta_writer:
    def __init__(self, fasta_file: str, cluster_file: str, output_path: str, bin_names: str, file_predicate: Callable[[str], bool]) -> None:
        self.fasta_file = fasta_file
        self.cluster_file = cluster_file
        self.file_predicate = file_predicate
        self.bin_names = bin_names
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
    def write_fasta(self, values: Dict[str, str]):
        n = 80 #how many chars per line
        print('writing new fasta file...')
        
        with open(self.output_path, 'x') as f:
            for binname, dna_string in tqdm(values.items()):
                f.write('>' + binname + '\n')
                for dna_line in [dna_string[i:i+n] for i in range(0, len(dna_string), n)]:
                    f.write(dna_line + '\n')
        
        
    def work(self, showClusters: bool = True):
        
        print(join(getcwd(), self.output_path))
        clusters = self.read_clusters()
        contigs = self.read_fasta()
        result = {}
    
        print("combining clusters...")
        for binname, contignames in tqdm(clusters.items()):
            result[self.bin_names+binname] = ''.join([contigs[name] for name in contignames])
        
        self.write_fasta(result)
        if showClusters:
            print(">Clusters:")
            for x in clusters.items():
                print(x)                
                
if __name__ == '__main__':
    print(sys.argv)
    writer = fasta_writer(sys.argv[1], sys.argv[2], sys.argv[3], \
        sys.argv[4] if len(sys.argv) >= 5 else 'bin_',\
        lambda x: x.startswith('cluster'))
    writer.work()