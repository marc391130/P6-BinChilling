from __future__ import annotations
import multiprocessing as mp
from multiprocessing import shared_memory
from multiprocessing.shared_memory import SharedMemory
from typing import Iterator, List, Tuple
import ctypes as c
from ClusterDomain import Cluster
import sys
import numpy as np


def mod(a: int, b: int) -> int:
    if b <= 0: raise Exception('Cannot mod with base less then or equal to 0')
    r = a % b
    return r + b if r < 0 else r

class HashTableIterator(Iterator):
    def __init__(self, hastable: SharedHashTable) -> None:
        self._table = hastable
        self._index = 0
        
    def __iter__(self) -> Iterator[Tuple[int, float]]:
        return HashTableIterator(self._table)
    
    def __next__(self) -> Tuple[int, float]:
        value = None
        while self._index < self._table.size:
            value = self._table.read_at_index(self._index)
            self._index += 1
            if value is not None:
                return value

        raise StopIteration()
                


class SharedHashTable:
    def __init__(self, size: int) -> None:
        self.size = size
        self.keys: mp.RawArray   = mp.RawArray(c.c_long, size)
        self.values: mp.RawArray = mp.RawArray(c.c_double, size)
    # def __init__(self, size: int, keys: np.ndarray, values: np.ndarray) -> None:
    #     self.size = size
    #     self.keys = keys
    #     self.values = values        
    
    def insert(self, key: int, value: float):
        key = 1 if key == 0 else key
        
        for hashf in [self.__hash_func_1__, self.__hash_func_2__, self.__hash_func_3__]:
            if self.__try_place__(hashf(key), key, value): return
        
        index = self.__hash_func_1__(key)
        for i in range(1, self.size):
            c_index = (index + i) % self.size
            if self.__try_place__(c_index, key, value): return
        
        raise Exception('Hashtable is full')
    
    def read_at_index(self, index: int) -> Tuple[int, float] or None:
        key, value = self.keys[index], self.values[index]
        if key == 0:
            return None
        return (key, value)
    
    def __try_place__(self, index: int, key: int, value: float) -> bool:
        key_v = self.keys[index]
        if key_v == 0:
            self.keys[index] = key
            self.values[index] = value
            return True
        if int(key_v) == key:
            raise KeyError('Key already exists')
        return False
    
    def remove(self, key: int) -> bool:
        for hashf in [self.__hash_func_1__, self.__hash_func_2__, self.__hash_func_3__]:
            index = hashf(key)
            if self.keys[index] == key:
                self.keys[index]   = 0
                self.values[index] = 0
                return True
        
        index = self.__hash_func_1__(key)
        for i in range(1, self.size):
            c_index = (index + i) % self.size
            if self.keys[c_index] == key:
                self.keys[index]   = 0
                self.values[index] = 0
                return True
            elif self.keys[c_index] == 0:
                return False
        
        return False
    
    def get(self, key: int, default_value: float) -> float:
        key = 1 if key == 0 else key
        for hashf in [self.__hash_func_1__, self.__hash_func_2__, self.__hash_func_3__]:
            index = hashf(key)
            result = self.__find_index__(index, key)
            if result is not None:
                return result
        #end for
        
        index = self.__hash_func_1__(key)
        result = self.__find_search__(key)
        return result if result is not None else default_value
    
    def __find_index__(self, index: int, key: int) -> float or None:
        key_v = self.keys[index]
        if key_v == 0:
            return None
        if key_v == key:
            return float(self.values[index])
        return None
    
    def __find_search__(self, key: int) -> float or None:
        has_seen_none = False 
        start_index = self.__hash_func_1__(key)
        for i in range(1, self.size):
            c_index = (start_index + i) % self.size
            key_v = self.keys[c_index]
            if key_v == 0:
                if has_seen_none:
                    return None
                has_seen_none = True
                continue
            elif key_v == key:
                return float(self.values[c_index])
        return None
    
    def test(self):
        pass
        # count = sum((1 for i in range(self.size) if self.keys[i] != 0))
        # print(f'{count} of {self.size} entries have value')
        
    def __hash_func_1__(self, key: int) -> int:
        return mod(key, self.size)
    
    def __hash_func_2__(self, key: int) -> int:
        n = c.c_ulong(key)
        ueven_const = 0x5555555555555555 #alternating byte e.g. 0101010101
        const = 17316035218449499591 #random uneven number
        def xorShift(n: int, i: int) -> int:
            return (n^(n>>i))
        
        value = const*xorShift(ueven_const*xorShift(n.value, 32), 32)
        return mod(value, self.size)
    
    def __hash_func_3__(self, key: int) -> int:
        x = key
        magic_number = c.c_ulonglong(0xbf58476d1ce4e5b9)
        x = (x ^ (x >> 30)) * magic_number.value
        x = (x ^ (x >> 27)) * magic_number.value
        x = x ^ (x >> 31)
        return mod(x, self.size)
    
    def __sizeof__(self):
        return sys.getsizeof(self.keys) + sys.getsizeof(self.values)
    

class UnevenNestedList:
    def __init__(self, starts: List[int], store: np.ndarray, copy_on_read = True) -> None:
        self._starts = starts
        self._store = store
        self.copy_on_read = copy_on_read
        #endfor
        
    @staticmethod
    def init_ram(cluster_lst: List[Cluster], total_elements: int = None) -> UnevenNestedList:
        count = len(cluster_lst)
        if total_elements is None:
            total_elements = sum((len(c) for c in cluster_lst))
            
        starts = np.empty(shape=count, dtype=np.intc)
        store = np.empty(shape=total_elements, dtype=np.longlong)
        
        index = 0
        for i, cluster in enumerate(cluster_lst):
            starts[i] = index
            for item in cluster:
                store[index] = hash(item)
                index += 1
        return UnevenNestedList(starts, store, copy_on_read=False)
    
    @staticmethod
    def init_combined(cluster_lst: List[Cluster], cluster_lst_2: List[Cluster]) -> Tuple[UnevenNestedList, int]:
        count = len(cluster_lst) + len(cluster_lst_2)
        
        total_elements = sum((len(c) for c in cluster_lst)) + sum((len(c) for c in cluster_lst_2))
            
        starts = np.empty(shape=count, dtype=np.intc)
        store = np.empty(shape=total_elements, dtype=np.longlong)
        
        index = 0
        max_i = 0
        for i, cluster in enumerate(cluster_lst):
            max_i = max(i, max_i)
            starts[i] = index
            for item in cluster:
                store[index] = hash(item)
                index += 1
        
        for i, cluster in enumerate(cluster_lst_2, start=len(cluster_lst)):
            starts[i] = index
            for item in cluster:
                store[index] = hash(item)
                index += 1
                
        return (UnevenNestedList(starts, store, copy_on_read=False), len(cluster_lst))
        
    @staticmethod
    def init_shared_memory(shared_mem_name: str, cluster_lst: List[Cluster], total_elements: int = None) -> UnevenNestedList:
        count = len(cluster_lst)
        if total_elements is None:
            total_elements = sum((len(c) for c in cluster_lst))
            
        starts = np.empty(shape=count, dtype=np.intc)
        shared_mem = SharedMemory(shared_mem_name, create=True, size=total_elements*sys.getsizeof(np.longlong))
        store = np.ndarray(shape=total_elements, dtype=np.longlong, buffer=shared_mem)
        
        index = 0
        for i, cluster in enumerate(cluster_lst):
            starts[i] = index
            for item in cluster:
                store[index] = hash(item)
                index += 1
        return UnevenNestedList(starts, store, copy_on_read=False)
        
    @staticmethod
    def init_memmap(filename:str , cluster_lst: List[Cluster], total_elements: int = None) -> UnevenNestedList:
        count = len(cluster_lst)
        if total_elements is None:
            total_elements = sum((len(c) for c in cluster_lst))
            
        _starts = np.empty(shape=count, dtype=np.intc)
        store_temp = np.empty(shape=total_elements, dtype=np.longlong)
        
        index = 0
        for i, cluster in enumerate(cluster_lst):
            _starts[i] = index
            for item in cluster:
                store_temp[index] = hash(item)
                index += 1
        
        store = np.memmap(filename=filename, dtype=np.longlong, shape=total_elements, mode='w+')
        store[:] = store_temp[:]
        store.flush()
        return UnevenNestedList(_starts, store)
                
    def __getitem__(self, cluster_index: int, item_index: int = None) -> int or List[int]:
        cluster = self.get_cluster(cluster_index)
        if item_index is None:
            return cluster
        return cluster[item_index]
    
    def get_item(self, cluster_index: int, item_index: int) -> int:
        cluster = self.get_cluster(cluster_index)
        return cluster[item_index]
    
    def get_cluster(self, index: int) -> List[int]:
        start_index = self._starts[index]
        end_index = self._starts[index+1] if (index < len(self._starts) -1) else (len(self._store))
        return np.array(self._store[start_index:end_index], dtype=np.longlong, copy=self.copy_on_read)
    
    def __len__(self) -> int:
        return len(self._starts)
    
    def __sizeof__(self):
        return sys.getsizeof(self._starts) + len(self._starts)*sys.getsizeof(c.c_int32) + sys.getsizeof(self._store) + len(self._store)*sys.getsizeof(c.c_longlong)