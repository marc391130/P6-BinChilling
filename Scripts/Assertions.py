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
