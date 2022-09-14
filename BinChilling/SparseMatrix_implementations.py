from __future__ import annotations
from email.policy import default
import itertools
from sqlite3 import NotSupportedError
from tqdm import tqdm
from typing import Iterable, Iterator, MutableMapping, Set, Tuple, Dict, Callable, TypeVar, Generic, List
import Assertions as Assert 
from math import sqrt

from SharedDatastructures_implementations import SharedHashTable

TK = TypeVar("TK")
TK2 = TypeVar("TK2")
TV = TypeVar("TV")


def SortKeysByHash(key1: TK, key2: TK) -> Tuple[TK, TK]:
    return (key1, key2) if key1.__hash__() <= key2.__hash__() else (key2, key1)


class HashIterator(Iterator[Tuple[TK, TV]]):
    def __init__(self, source: Dict[Tuple[TK, TK], TV], keysort : Callable[[TK, TK], Tuple[TK, TK]],\
        pivot: TK, control: Iterator[TK], second_iterator: Iterator[TK] = None) -> None:
        self.source = source
        self.pivot = pivot
        self.control = control if second_iterator is None else itertools.chain(control, second_iterator)
        self.keysort = keysort
        
    def __iter__(self) -> Iterator[Tuple[TK, TV]]:
        return self
        
    def __next__(self) -> Tuple[TK, TV]:
        tupkey = self.keysort(self.pivot, self.control.__next__())
        return self.source[tupkey]

class SparseDictHashMatrix(MutableMapping[Tuple[TK, TK2], TV]):
    def __init__(self, keysort: Callable[[TK, TK2], Tuple[TK, TK2]] = None, default_value = None,\
        sparse_value_predicate: Callable[[TV], bool] = None) -> None:
        #sparses the value if the predicate returns true 
        
        self.__internal__ : Dict[TK, Dict[TK2, TV]] = dict()
        self.keysort = keysort if keysort is not None else lambda x,y: (x,y)
        self.__default__ = default_value
        self.__sparse_filter__ = sparse_value_predicate if sparse_value_predicate is not None\
            else lambda x: x == self.__default__
        
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK2) -> TV:
        k1, k2 = self.keysort(key1, key2)
        return self.get( (k1, k2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK2]]:
        return self.__internal__.keys()
    
    def values(self) -> Iterable[TV]:
        return self.__internal__.values()
    
    def items(self) -> Iterable[Tuple[TK, Dict[TK2, TV]]]:
        return self.__internal__.items()
    
    def get(self, __k: Tuple[TK, TK2]) -> TV:
        k1, k2 = self.keysort(__k[0], __k[1])
        return self.__internal__.get(k1, {}).get(k2, self.__default__)
    
    def __getitem__(self, __k: Tuple[TK, TK2]) -> TV:
        return self.get(__k)
    
    def get_row(self, __k: TK) -> Dict[TK2, TV]:
        return self.__internal__.get(__k, {})
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK2], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK2], __v: TV) -> None:
        if self.__sparse_filter__(__v): return
        
        k1, k2 = self.keysort(__k[0], __k[1])
        if k1 not in self.__internal__: self.__internal__[k1] = {}
        self.__internal__[k1][k2] = __v
    
    def set_entry(self, k1: TK, k2: TK2, v: TV) -> None:
        return self.set( (k1, k2), v )
    
    def set_dict(self, key: TK, dct: Dict[TK2, TV]) -> None:
        for other_key, value in dct.items():
            self.set( self.keysort( key, other_key ), value)

    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK2]) -> TV:
        k1, k2 = self.keysort(__v[0], __v[1])
        return self.__internal__[k1].pop(k2)
    
    def pop(self, __k: Tuple[TK, TK2]) -> None:
        return self.__delitem__(__k)
    
    def pop_entry(self, key1: TK, key2: TK2) -> None:
        return self.__delitem__( self.keysort(key1, key2) )

    def pop_contain(self, key: TK) -> None: #will always return none and never throw
        for key, inner_dct in self.__internal__:
            inner_dct.pop(key, None)
        self.__internal__.pop(key, None)
        
    def pop_set(self, remove_keys: Set[TK or TK2]) -> None:
        
        if len(remove_keys) == 0: return
        pk_remove_set = []
        for pk, column in self.__internal__.items():
            if pk in remove_keys: pk_remove_set.append(pk)
            else:
                intersection = remove_keys.intersection(column.keys())
                for sk in intersection:
                    column.pop(sk)
        
        for pk_to_remove in pk_remove_set:
            self.__internal__.pop(pk_to_remove)
        
    
    #UTILITY FUNCTIONS
    def has_row_key(self, key: TK) -> bool:
        return key in self.__internal__
    
    def has_tuple(self, key: Tuple[TK, TK2]) -> bool:
        k1, k2 = self.keysort(key[0], key[1])
        return self.has_entry(k1, k2)
    
    def has_entry(self, k1: TK, k2: TK2) -> bool:
        k1, k2 = self.keysort(k1, k2)
        return k1 in self.__internal__ and k2 in self.__internal__[k1]
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        return False
    
    def count_entries(self) -> int:
        return sum( (len(row_dct) for _, row_dct in self.__internal__.items()) )
    
    def __len__(self) -> int:
        return len(self.__internal__)
    
    def __iter__(self) -> Iterator[Tuple[TK, Dict[TK2, TV]]]:
        return self.__internal__.__iter__()
        
class DoubleSparseDictHashMatrix(MutableMapping[Tuple[TK, TK2], TV]):
    def __init__(self, default_value = None, sparse_value_predicate: Callable[[TV], bool] = None) -> None:
        #sparses the value if the predicate returns true 
        
        self.__internal_row__ : Dict[TK, Dict[TK2, TV]] = dict()
        self.__internal_column__ : Dict[TK2, Dict[TK, TV]] = dict()
        self.__default__ = default_value
        self.__sparse_filter__ = sparse_value_predicate if sparse_value_predicate is not None\
            else lambda x: x == self.__default__
        
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK2) -> TV:
        return self.get( (key1, key2) )
    
    
    def items(self) -> Iterable[Tuple[TK, Dict[TK, TV]]]:
        return self.__internal_row__.items()
    
    def get(self, __k: Tuple[TK, TK2]) -> TV:
        k1, k2 = __k
        return self.__internal_row__.get(k1, {}).get(k2, self.__default__)
    
    def __getitem__(self, __k: Tuple[TK, TK2]) -> TV:
        return self.get(__k)
    
    def get_row(self, __k: TK) -> Dict[TK2, TV]:
        return dict(self.__internal_row__.get(__k, {}))
    
    def get_column(self, __k: TK2) -> Dict[TK, TV] or None:
        return dict(self.__internal_column__.get(__k, {}))
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK2], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK2], __v: TV) -> None:
        if self.__sparse_filter__(__v): return
        k1, k2 = __k
        if k1 not in self.__internal_row__: self.__internal_row__[k1] = {}
        self.__internal_row__[k1][k2] = __v
        if k2 not in self.__internal_column__: self.__internal_column__[k2] = {}
        self.__internal_column__[k2][k1] = __v
    
    def set_row(self, key: TK, dct: Dict[TK2, TV]) -> None:
        for other_key, value in dct.items():
            self.set( (key, other_key) , value)

    def set_column(self, key: TK2, dct: Dict[TK, TV]) -> None:
        for other_key, value in dct.items():
            self.set( ( other_key, key), value)
    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK2]) -> TV:
        k1, k2 = __v
        x1 = self.__internal_row__[k1].pop(k2, self.__default__)
        x2 = self.__internal_column__[k2].pop(k1, self.__default__)
        return x1
    
    def pop(self, __k: Tuple[TK, TK2]) -> None:
        return self.__delitem__(__k)
    
    def pop_entry(self, key1: TK, key2: TK2) -> None:
        return self.__delitem__( (key1, key2) )
    
    #UTILITY FUNCTIONS
    def __len__(self) -> Tuple[int, int]:
        return (len(self.__internal_row__), len(self.__internal_column__))
    
    def has_row_key(self, key: TK) -> bool:
        return key in self.__internal_row__
    
    def has_column_key(self, key: TK2) -> bool:
        return key in self.__internal_column__
    
    def has_tuple(self, key: Tuple[TK, TK2]) -> bool:
        k1, k2 = key
        return self.has_entry(k1, k2)
    
    def has_entry(self, k1: TK, k2: TK2) -> bool:
        return k1 in self.__internal_row__ and k2 in self.__internal_row__.get(k1, {})
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        return False
    
    def __iter__(self) -> Iterator[Tuple[TK, Dict[TK2, TV]]]:
        return self.__internal_row__.__iter__()
    
    
class SparseTupleHashMatrix(MutableMapping[Tuple[TK, TK], TV]):
    def __init__(self, keysort : Callable[[TK, TK], Tuple[TK, TK]] = None, default_value = 0.0) -> None:
        self.__internal__ : Dict[Tuple[TK, TK], TV] = dict()
        self.default_value = default_value
        self.keysort = keysort if keysort is not None else lambda x, y: (x, y) 
    
    def sortTuple(self, tup: Tuple[TK, TK]) -> Tuple[TK, TK]:
        return self.keysort(tup[0], tup[1])
    
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK) -> TV:
        return self.get( self.keysort(key1, key2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK]]:
        return self.__internal__.keys()
    
    def values(self) -> Iterable[TV]:
        return self.__internal__.values()
    
    def items(self) -> Iterable[Tuple[TK, TK], TV]:
        return self.__internal__.items()
    
    def get(self, __k: Tuple[TK, TK]) -> TV:
        return self.__internal__.get(self.sortTuple(__k), self.default_value)
    
    def __getitem__(self, __k: Tuple[TK, TK]) -> TV:
        return self.get(__k)
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK], __v: TV) -> None:
        tup = self.sortTuple(__k)
        self.__internal__[tup] = __v
    
    def set_dict(self, key: TK, dct: Dict[TK, TV] ):
        for other_key, value in dct.items():
            self.set( (key, other_key), value )
    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK]) -> None:
        return self.pop(__v)
    
    def pop(self, __k: Tuple[TK, TK]) -> None:
        self.__internal__.pop(__k)
    
    def pop_entry(self, key1: TK, key2: TK) -> None:
        tup = self.keysort(key1, key2)
        Assert.assert_key_exists(tup, self.__internal__)
        self.__internal__.pop(tup)

    def pop_contain(self, key: TK) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] is key or tup[1] is key: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    def pop_set(self, keys: set[TK]) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] in keys or tup[1] in keys: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    #UTILITY FUNCTIONS
    def has_entry(self, key: Tuple[TK, TK]) -> bool:
        tup = self.sortTuple(key)
        return tup in self.__internal__
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        return False
    
    def __len__(self) -> int:
        return len(self.__internal__)
    
    def __iter__(self) -> Iterator[Tuple[TK, TK], TV]:
        return self.__internal__.__iter__()
    
class SuperSparseTupleHashMatrix(MutableMapping[Tuple[TK, TK], TV]):
    def __init__(self, keysort : Callable[[TK, TK], Tuple[TK, TK]] = None, default_value = 0.0) -> None:
        self.__internal__ : Dict[int, TV] = dict()
        self.default_value = default_value
        self.keysort = keysort if keysort is not None else lambda x, y: (x, y) 
    
    def sortTuple(self, tup: Tuple[TK, TK]) -> Tuple[TK, TK]:
        return self.keysort(tup[0], tup[1])
    
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK) -> TV:
        return self.get( self.keysort(key1, key2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK]]:
        raise NotSupportedError(".keys not supported on SuperSparseTupleHashMatrix")
    
    def values(self) -> Iterable[TV]:
        return self.__internal__.values()
    
    
    def get(self, __k: Tuple[TK, TK]) -> TV:
        index = hash(self.sortTuple(__k))
        return self.__internal__.get(index, self.default_value)
    
    def __getitem__(self, __k: Tuple[TK, TK]) -> TV:
        return self.get(__k)
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK], __v: TV) -> None:
        index = hash(self.sortTuple(__k))
        self.__internal__[index] = __v
    
    def set_dict(self, key: TK, dct: Dict[TK, TV] ):
        for other_key, value in dct.items():
            self.set( (key, other_key), value )
    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK]) -> None:
        return self.pop(__v)
    
    def pop(self, __k: Tuple[TK, TK]) -> None:
        self.__internal__.pop(__k)
    
    def pop_entry(self, key1: TK, key2: TK) -> None:
        tup = self.keysort(key1, key2)
        Assert.assert_key_exists(tup, self.__internal__)
        self.__internal__.pop(tup)

    def pop_contain(self, key: TK) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] is key or tup[1] is key: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    def pop_set(self, keys: set[TK]) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] in keys or tup[1] in keys: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    #UTILITY FUNCTIONS
    def has_entry(self, key: Tuple[TK, TK]) -> bool:
        index = hash(self.sortTuple(key))
        return index in self.__internal__
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        if type(__o) is int:
            return __o in self.__internal__
        return False
    
    def __len__(self) -> int:
        return len(self.__internal__)
    
    def __iter__(self) -> Iterator[Tuple[TK, TK], TV]:
        return self.__internal__.__iter__()


class SparseTupleHashMatrix(MutableMapping[Tuple[TK, TK], TV]):
    def __init__(self, keysort : Callable[[TK, TK], Tuple[TK, TK]] = None, default_value = 0.0) -> None:
        self.__internal__ : Dict[Tuple[TK, TK], TV] = dict()
        self.default_value = default_value
        self.keysort = keysort if keysort is not None else lambda x, y: (x, y) 
    
    def sortTuple(self, tup: Tuple[TK, TK]) -> Tuple[TK, TK]:
        return self.keysort(tup[0], tup[1])
    
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK) -> TV:
        return self.get( self.keysort(key1, key2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK]]:
        return self.__internal__.keys()
    
    def values(self) -> Iterable[TV]:
        return self.__internal__.values()
    
    def items(self) -> Iterable[Tuple[TK, TK], TV]:
        return self.__internal__.items()
    
    def get(self, __k: Tuple[TK, TK]) -> TV:
        return self.__internal__.get(self.sortTuple(__k), self.default_value)
    
    def __getitem__(self, __k: Tuple[TK, TK]) -> TV:
        return self.get(__k)
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK], __v: TV) -> None:
        tup = self.sortTuple(__k)
        self.__internal__[tup] = __v
    
    def set_dict(self, key: TK, dct: Dict[TK, TV] ):
        for other_key, value in dct.items():
            self.set( (key, other_key), value )
    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK]) -> None:
        return self.pop(__v)
    
    def pop(self, __k: Tuple[TK, TK]) -> None:
        self.__internal__.pop(__k)
    
    def pop_entry(self, key1: TK, key2: TK) -> None:
        tup = self.keysort(key1, key2)
        Assert.assert_key_exists(tup, self.__internal__)
        self.__internal__.pop(tup)

    def pop_contain(self, key: TK) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] is key or tup[1] is key: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    def pop_set(self, keys: set[TK]) -> None:
        to_remove = set()
        for tup in self.__internal__.keys():
            if tup[0] in keys or tup[1] in keys: to_remove.add(tup)
        for el in to_remove:
            self.__internal__.pop(el)
    
    #UTILITY FUNCTIONS
    def has_entry(self, key: Tuple[TK, TK]) -> bool:
        tup = self.sortTuple(key)
        return tup in self.__internal__
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        return False
    
    def __len__(self) -> int:
        return len(self.__internal__)
    
    def __iter__(self) -> Iterator[Tuple[TK, TK], TV]:
        return self.__internal__.__iter__()
    
class HashedMatrix(MutableMapping[Tuple[TK, TK], TV]):
    def __init__(self, key_hash: Callable[[object], int], base: SharedHashTable) -> None:
        self.base = base
        self.key_hash = key_hash 
    
    #READ FUNCTIONS
    def getEntry(self, key1: TK, key2: TK) -> TV:
        return self.get( (key1, key2) )
    
    def keys(self) -> Iterable[Tuple[TK, TK]]:
        raise NotSupportedError(".keys not supported on HashedMatrix")
    
    def values(self) -> Iterable[TV]:
        raise NotSupportedError(".values not supported on HashedMatrix")
    
    
    def get(self, __k: Tuple[TK, TK], default_value: TV = None) -> TV:
        hash_value = self.key_hash( __k )
        return self.base.get(hash_value, default_value)
    
    def __getitem__(self, __k: Tuple[TK, TK]) -> TV:
        return self.get(__k)
    
    #SET FUNCTIONS
    def __setitem__(self, __k: Tuple[TK, TK], __v: TV) -> None:
        self.set(__k, __v)
    
    def set(self, __k: Tuple[TK, TK], __v: TV) -> None:
        index = self.key_hash(__k)
        self.base.insert(index, __v)
    
    def set_dict(self, key: TK, dct: Dict[TK, TV] ):
        for other_key, value in dct.items():
            self.set( (key, other_key), value )
    
    #DELETE FuNCTIONS
    def __delitem__(self, __v: Tuple[TK, TK]) -> None:
        return self.pop(__v)
    
    def pop(self, __k: Tuple[TK, TK]) -> None:
        index = self.key_hash(__k)
        self.base.remove(index)
    
    def pop_entry(self, key1: TK, key2: TK) -> None:
        self.pop( (key1, key2) )

    #UTILITY FUNCTIONS
    def has_entry(self, key: Tuple[TK, TK]) -> bool:
        result = self.get(key, default_value=None)
        return result is not None
        
    
    def __contains__(self, __o: object) -> bool:
        if type(__o) is Tuple:
            return self.has_entry(__o)
        if type(__o) is int:
            return self.base.get(__o, default_value=None) is not None
        return False
    
    def count_entries(self) -> int:
        return sum((1 if x != 0 else 0 for x in self.base.keys))
    
    def __len__(self) -> int:
        return self.base.size
    
    def __iter__(self) -> Iterator[Tuple[TK, TK], TV]:
        return self.__internal__.__iter__()