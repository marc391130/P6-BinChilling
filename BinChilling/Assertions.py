from typing import Callable, Dict, List, Collection

def assert_fail(msg: str = None):
    err = "An assertion failed" if msg is None else msg 
    raise AssertionError(err)

def assert_not_none(item: object):
    if item is None:
        raise AssertionError(f"item {item} is None")

def assert_index_exists(idx: int, lst: List) -> None:
    if idx not in lst:
        raise AssertionError(f"Index {idx} not in {str(lst)}")

def assert_item_in_list(lst: List, item: any):
    if item not in lst:
        raise AssertionError(f"Item {item} not in list")

def assert_item_not_in_collection(lst: Collection, item: any):
    if item in lst:
        raise AssertionError(f"Item '{item}' in collection")

def assert_key_exists(key, dct: Dict) -> None:
    if key not in dct:
        raise AssertionError(f"Key {key} not in {str(dct)}")
    
def assert_key_not_exists(key, dct: Dict) -> None:
    if key in dct:
        raise AssertionError(f"Key {key} in {str(dct)}")

def assert_partition_content(data: List, item) -> None:
    if item not in data:
        raise AssertionError(f"Item {item} not in data")

def assert_unique(collection: Collection):
    if len(collection) <= 1:
        return
    
    
    if len(collection) > len(set(collection)):
        raise AssertionError("Collection is not consisting of unique elements")

def assert_partion_set_content(partition_set: List[Dict]) -> None:
    raf = None
    for partition in partition_set:
        if raf is None:
            raf = partition.__data__

        if partition.__data__ is not raf:
            raise AssertionError("Data is not the same in partitionset and partition")

def assert_min_list_len(lst: List, min_len: int):
    if len(lst) < min_len:
        raise AssertionError(f"list has length {len(lst)}, but should have a minimum of {min_len}")

def assert_list_nonempty(lst: List) -> None:
    if len(lst) <= 0:
        raise AssertionError("list is empty, when it shouldnt be")

def assert_new_cluster(object1, object2) -> None:
    if object1 is None and object2 is not None:
        raise AssertionError("Object1 is none, while object2 is not")
    elif object1 is not None and object2 is None:
        raise AssertionError("Object2 is none, while object1 is not")
    elif object1 is not None and object1 is object2:
        raise AssertionError("Trying to make a new cluster, with two of the same references, not allowed!")

def assert_objects_of_list_is_type(lst, target_type) -> None:
    for item in lst: 
        if type(item) != target_type:
            raise AssertionError(f"Object {item} is not of correct type: {target_type}, but is instead of type: {type(item)}")

def assert_only_one_of_each_element_in_cluster_list(elements, cluster_lst) -> None:
    for element in elements:
        elements_found = 0
        for cluster in cluster_lst:
            if element in cluster:
                elements_found += 1
    
        if elements_found > 1:
            raise AssertionError("Too many instances of element:", element)

def assert_all_elements_are_in_cluster_lst(elements, cluster_lst) -> None:
    item_set = set()

    for cluster in cluster_lst:
        for item in cluster:
            if item not in item_set: 
                item_set.add(item)
    
    for element in elements:
        if element not in item_set:
            raise AssertionError("Missing element in item_set")

def assert_lambda(assertion_predicate: Callable[[], bool], error_msg:str =None):
    if assertion_predicate() is False:
        msg = f"'{assertion_predicate}' assertion failed" if error_msg is None else error_msg
        raise AssertionError(msg)
    

def assert_equal(value1, value2) -> None:
    if value1 != value2:
        raise AssertionError(f"{value1} != {value2}")

def assert_not_equal(value1, value2) -> None:
    if value1 == value2:
        raise AssertionError(f"{value1} == {value2}")
    
def assert_in_range(value, minv, maxv) -> None:
    minv, maxv = (minv, maxv) if minv <= maxv else (maxv, minv)
    
    if minv <= value <= maxv is False:
        raise AssertionError(f"'{value}' is not in range {minv} - {maxv}")
        