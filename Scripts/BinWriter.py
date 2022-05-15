from os import listdir
from os.path import join
from typing import Callable, List, Dict, Tuple
import sys
from tqdm import tqdm

ATCG = set(['A', 'T', 'C', 'G', '\n'])

class bin_writer:
    def __init__(self, fasta_file: str, output_path: str, cluster_file: str or None) -> None:
        self.fasta_file = fasta_file
        self.cluster_file = cluster_file
        self.output_path = output_path
        
    # key is contig name, value is contig string 
    def read_fasta(self) -> Tuple[Dict[str, List[str]], Dict[str, bool]]:
        current = ''
        result: Dict[str, List[str]] = {}
        is_ATCG : Dict[str, bool] = {}
        
        def contig_cleaner(contingname: str) -> str:
            return contingname.replace('>', '').replace('\n', '')
        
        print("reading fasta...")
        
        with open(self.fasta_file, 'r') as f:
            for line in tqdm(f.readlines()):
                line: str
                if line.startswith('>'):
                    current = contig_cleaner(line)
                    result[current] = []
                    is_ATCG[current] = True
                else:
                    result[current].append(line)
                    if all(char in ATCG for char in line) is False:
                        is_ATCG[current] = False
                    
        # print('Cleaning dna string for escape characters')
        # for bin, dna in result.items():
        #     result[bin] = dna.replace('\n', '')
        return result, is_ATCG
    
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
        outputpath = outputpath if outputpath.endswith(".fasta") else outputpath + ".fasta"
        
        with open(outputpath, 'x') as f:
            for binname, dna_string in values.items():
                f.write('>' + binname + '\n')
                for dna_line in [dna_string[i:i+n] for i in range(0, len(dna_string), n)]:
                    f.write(dna_line + '\n')
    
    def write_fasta_fast(self, outputpath: str, values: Dict[str, List[str]]):
        with open(outputpath, 'x') as f:
            for binname, dna_string_lst in values.items():
                f.write('>' + binname + '\n')
                for dna_line in dna_string_lst:
                    f.write(dna_line)

    def write_bins(self, values: Dict[str, Dict[str, List[str]]], keep_lst: Dict[str, bool] = None, extension: str = 'fasta'):
        print("writing bins...")
        keep_lst = keep_lst if keep_lst is not None else {}
        count = 0
        for clustername, cluster in tqdm(values.items()):
            if clustername in keep_lst and keep_lst[clustername] is False: 
                count += 1
                print(f'skipping {clustername}. ({count} / {len(keep_lst)})')
                continue
            name = str(clustername) + '.' + extension
            self.write_fasta_fast(join(self.output_path, name), cluster)
        
    def read_cluster_from_fasta(self) -> Dict[str, List[str]]:
        result = {}
        print("reading clusters...")
        def contig_cleaner(contingname: str) -> str:
            return contingname.replace('>', '').replace('\n', '')
        
        with open(self.fasta_file, 'r') as f:
            for line in tqdm(f.readlines()):
                line: str
                if line.startswith('>'):
                    clean = contig_cleaner(line)
                    result[clean] = [clean]
        return result  

    
    
    def work(self, showClusters: bool = False):
        
        clusters = self.read_clusters() if self.cluster_file is not None else self.read_cluster_from_fasta()
        contigs, is_ATCG = self.read_fasta()
        
        #split contigs into clusters
        print("split contigs into clusters...")
        result: Dict[str, Dict[str, List[str]]] = {binname: {name: contigs[name] for name in contignames} \
            for binname, contignames in tqdm(clusters.items())} #cluster name to dict of contigname to contingstring 
        
        self.write_bins(result, is_ATCG, extension=self.fasta_file.split('.')[-1])
        
        if showClusters:
            for k, v in result.items():
                print(f">{k}") 
            for k2, v2 in v.items():
                print(f"{k}_{k2} len: {len(v2)}")               
                
if __name__ == '__main__':
    print(sys.argv)
    if len(sys.argv) == 4:
        print('Fasta_file, output_path, partition_file')
    elif len(sys.argv) == 3:
        print('Fasta_file, output_path')
    else:
        print('Fasta_file, output_path, partition_file (optional)')
        sys.exit()
    
    writer = bin_writer(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) >= 4 else None)
    writer.work()