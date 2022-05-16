import sys
import os
from tqdm import tqdm

def main():
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    
    count, max_count = 1, len(os.listdir(input_folder))
    
    print('Running mover')
    for file in os.listdir(input_folder):
        filepath = os.path.join(input_folder, file)
        out_folder = os.path.join(output_folder, file.split('.')[0])
        os.mkdir(out_folder)
        os.rename(filepath, os.path.join(out_folder, file))
        
    

if __name__ == '__main__':
    main()