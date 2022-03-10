import sys
from typing import Dict

filepath = sys.argv[1]

def read_file(filepath: str) -> Dict[str, int]:
    answer_list = {'P': 0, 'L': 0, 'S': 0}
    P_dict = {}

    with open(filepath, 'r') as file:
        for line in file.readlines():
            if len(line) > 0:
                s = line[0]
                
                if s == 'P':
                    contig_name = line.split('\t')[1].split('_')[1]
                    if contig_name in P_dict:
                        P_dict[contig_name] += 1
                    else:
                        P_dict[contig_name] = 1

                if s in answer_list:
                    answer_list[s] += 1

    return P_dict

def find_missing(P_dict):
    for i in range(1150):
        if str(i) not in P_dict:
            print(i)
    
find_missing(read_file(filepath))


