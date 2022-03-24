from typing import Dict, List

def assert_index_exists(idx: int, lst: List) -> None:
    if idx not in lst:
        raise Exception(f"Index {idx} not in {str(lst)}")

def assert_item_in_list(lst: List, item: any):
    if item not in lst:
        raise Exception(f"Item {item} not in list")

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

def assert_list_has_atleast_len(lst: List, min_len: int):
    if len(lst) <= min_len:
        raise Exception(f"list has length {len(lst)}, but should have a minimum of {min_len}")

def assert_list_nonempty(lst: List) -> None:
    if len(lst) <= 0:
        raise Exception("list is empty, when it shouldnt be")

def assert_new_cluster(object1, object2) -> None:
    if object1 is None and object2 is not None:
        raise Exception("Object1 is none, while object2 is not")
    elif object1 is not None and object2 is None:
        raise Exception("Object2 is none, while object1 is not")
    elif object1 is not None and object1 is object2:
        raise Exception("Trying to make a new cluster, with two of the same references, not allowed!")

def assert_objects_of_list_is_type(lst, target_type) -> None:
    for item in lst: 
        if type(item) != target_type:
            raise Exception(f"Object {item} is not of correct type: {target_type}, but is instead of type: {type(item)}")

def assert_only_one_of_each_element_in_cluster_list(elements, cluster_lst) -> None:
    for element in elements:
        elements_found = 0
        for cluster in cluster_lst:
            if element in cluster:
                elements_found += 1
    
        if elements_found > 1:
            raise Exception("Too many instances of element:", element)

def assert_all_elements_are_in_cluster_lst(elements, cluster_lst) -> None:
    item_set = set()

    for cluster in cluster_lst:
        for item in cluster:
            if item not in item_set: 
                item_set.add(item)
    
    for element in elements:
        if element not in item_set:
            raise Exception("Missing element in item_set")
