import sys
import os
from tqdm import tqdm

def run_checkm(cpus: int, type: str, filepath: str, Bacteria_file: str, output_folder: str):
    os.system(f'checkm analyze -t {cpus} -x {type} {Bacteria_file} {filepath} {output_folder}')

def main():
    input_folder = [x[0] for x in os.walk(sys.argv[1])]
    output_folder = sys.argv[2]
    
    for folder in input_folder:
        pass
        
    

if __name__ == '__main__':
    main()