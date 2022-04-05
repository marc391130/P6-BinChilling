from itertools import groupby
from typing import List, Dict, Tuple
from tqdm import tqdm
import sys


def ClusterEquality(lst1: List, lst2: List) -> bool:
    if len(lst1) != len(lst2):
        return False
    
    for item in lst1:
        if item not in lst2:
            return False
    return True


class BinEqualityAsserter:
    def __init__(self, cluster_path1: str, cluster_path2: str) -> None:
        self.filepath1 = cluster_path1
        self.filepath2 = cluster_path2
        
    #EDGE, cluster
    def __read_file__(self, filepath: str) -> List[Tuple[str, str]]:
        result = []
        print(f"reading file { filepath } ...")
        with open(filepath, 'r') as f:
            for line in tqdm(f):
                split = line.replace('\n', '').split('\t')
                result.append((split[0], split[1]))
        return result
    
    def run(self) -> Tuple[bool, str]:
        group_func = lambda x: x[0]
        cluster1 = groupby(self.__read_file__(self.filepath1), group_func)
        cluster2 = groupby(self.__read_file__(self.filepath2), group_func)
        len1 = len(list(cluster1))
        len2 = len(list(cluster2))
        if len1 != len2:
            return False, f"Not same size. cluster1 contains {len1} elements, but cluster2 contains {len2}"
        
        print("checking equality...")
        for key1, group1 in tqdm(cluster1):
            if any(ClusterEquality(group1, group2) for key2, group2 in cluster2) is False:
                return False, f"Cluster { key1 } i første fil stemmer ikke overens med anden fil"
        
        return True, "Equal"

if __name__ == "__main__":
    print(sys.argv)
    asserter = BinEqualityAsserter(sys.argv[1], sys.argv[2])
    equal, msg = asserter.run()
    
    if equal:
        print("The files are equal")
    else:
        print("The files are not equal. \nfailed on: " + msg) 