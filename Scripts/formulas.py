from math import sqrt
from typing import List

def cluster_simularity(cluster1: List, cluster2: List, intersection_len: int,  total_elements: int) -> float:    
    counter = intersection_len - ((len(cluster1) * len(cluster2)) / total_elements)
    divisor = sqrt(len(cluster1) * len(cluster2) * (1 - (len(cluster1) / total_elements)) * (1 - (len(cluster2) / total_elements)))

    return counter / divisor if divisor != 0 else 0
