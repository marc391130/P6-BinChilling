from typing import Dict, List

def assert_index_exists(idx: int, lst: List) -> None:
    if idx not in lst:
        raise Exception(f"Index {idx} not in {str(lst)}")

def assert_key_exists(key, dct: Dict) -> None:
    if key not in dct:
        raise Exception(f"Key {key} not in {str(dct)}")

def assert_partition_content(data: List, item) -> None:
    if item not in data:
        raise Exception(f"Item {item} not in data")
        

def assert_partion_set_content(partition_set: List[Dict]) -> None:
    raf = None
    for partition in partition_set:
        if raf is None:
            raf = partition.__data__

        if partition.__data__ is not raf:
            raise Exception("Data is not the same in partitionset and partition")

def assert_new_cluster(object1, object2) -> None:
    if object1 is None and object2 is not None:
        raise Exception("Object1 is none, while object2 is not")
    elif object1 is not None and object2 is None:
        raise Exception("Object2 is none, while object1 is not")
    elif object1 is not None and object1 is object2:
        raise Exception("Trying to make a new cluster, with two of the same references, not allowed!")

