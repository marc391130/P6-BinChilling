from tqdm import tqdm
from cython.view cimport array
from cython.parallel cimport prange, parallel

cdef extern from "math.h":
    int floor(double) nogil

cdef int mod(int a, int b) nogil:
    cdef int r
    if b <= 0: 
        with gil:
            raise Exception('Cannot mod with base less then or equal to 0')
    r = a % b
    return r + b if r < 0 else r

cdef double[:] init_double_array(int size, double init_value):
    cdef int i
    cdef double[:] temp = array(shape=(size,), itemsize=sizeof(double), format='d')
    for i in range(size):
        temp[i] = init_value
    return temp

cdef long[:] init_long_array(int size, long init_value):
    cdef int i
    cdef long[:] temp = array(shape=(size,), itemsize=sizeof(long), format='l')
    for i in range(size):
        temp[i] = init_value
    return temp

cdef class HashTable:
    cdef int size
    cdef long[:] keys
    cdef double[:] values

    def __init__(self, size: int) -> None:
        self.size = size
        self.keys = init_long_array(size, 0)
        self.values = init_double_array(size, 0.0)
    # def __init__(self, size: int, keys: np.ndarray, values: np.ndarray) -> None:
    #     self.size = size
    #     self.keys = keys
    #     self.values = values        
    
    cdef insert(self, int key, double value):
        cdef int index
        key = 1 if key == 0 else key
        
        index = self.__hash_func_1__(key)
        if self.__try_place__(index, key, value): return
        if self.__try_place__(self.__hash_func_2__(key), key, value): return
        if self.__try_place__(self.__hash_func_3__(key), key, value): return

        for i in range(1, self.size):
            c_index = (index + i) % self.size
            if self.__try_place__(c_index, key, value): return
        
        raise Exception('Hashtable is full')
    
    cdef bint __try_place__(self, int index, long key, double value):
        cdef long key_v
        key_v = self.keys[index]
        if key_v == 0:
            self.keys[index] = key
            self.values[index] = value
            return True
        if key_v == key:
            raise KeyError('Key already exists')
        return False
    
    cdef double get(self, long key, double default_value) nogil:
        cdef double result
        cdef int index
        key = 1 if key == 0 else key

        index = self.__hash_func_1__(key)
        result = self.__find_index__(index, key)
        if result >= 0: return result


        index = self.__hash_func_2__(key)
        result = self.__find_index__(index, key)
        if result >= 0: return result

        index = self.__hash_func_3__(key)
        result = self.__find_index__(index, key)
        if result >= 0: return result
        
        index = self.__hash_func_1__(key)
        result = self.__find_search__(key)
        return result if result >= 0.0 else default_value
    
    cdef double __find_index__(self, int index, long key) nogil:
        cdef long key_v
        key_v = self.keys[index]
        if key_v == 0:
            return -1.0
        if key_v == key:
            return self.values[index]
        return -1.0
    
    cdef double __find_search__(self, long key) nogil:
        cdef bint has_seen_none = False 
        cdef int i, start_index, c_index
        cdef long key_v

        start_index = self.__hash_func_1__(key)
        for i in range(1, self.size):
            c_index = (start_index + i) % self.size
            key_v = self.keys[c_index]
            if key_v == 0:
                if has_seen_none:
                    return -1.0
                has_seen_none = True
                continue
            elif key_v == key:
                return self.values[c_index]
        return -1.0
    
        
    cdef int __hash_func_1__(self, long key) nogil:
        return mod(key, self.size)
    
    cdef int xorShift(self, int n,int i) nogil:
            return (n^(n>>i))

    cdef int __hash_func_2__(self, long key) nogil:
        cdef unsigned long rvalue, value
        cdef unsigned long n = <unsigned long> key 
        #alternating byte e.g. 0101010101
        cdef unsigned long ueven_const = 0x5555555555555555 
        rvalue = 17316035218449499591 #random uneven number
        
        value = rvalue*self.xorShift(ueven_const*self.xorShift(n, 32), 32)
        return mod(value, self.size)
    
    cdef int __hash_func_3__(self, long key) nogil:
        cdef long x
        cdef unsigned long long magic_number
        x = key
        magic_number = 0xbf58476d1ce4e5b9
        x = (x ^ (x >> 30)) * magic_number
        x = (x ^ (x >> 27)) * magic_number
        x = x ^ (x >> 31)
        return mod(x, self.size)
    


cdef HashTable co_matrix = None

cpdef int calc(int a, int b):
    return a + b


cdef long hash(long v1, long v2) nogil:
    cdef unsigned long mult = 0xf4243
    cdef unsigned long x = <unsigned long> 0x345678
    cdef int ln = 2
    cdef unsigned long y

    y = <unsigned long>v1
    x = (x ^ y) * mult
    mult += 82520UL + 4
    
    y = <unsigned long>v2
    x = (x ^ y) * mult

    x += 97531UL

    return <long>x

cpdef void init_coassosiation(dict values, int size):
    global co_matrix
    cdef int k1, k2, hash_value
    cdef tuple tup_key
    cdef double value
    co_matrix = HashTable(size)

    for tup_key, value in values.items():
        k1, k2 = tup_key
        hash_value = hash(k1, k2)
        co_matrix.insert(hash_value, value)
        print(str(hash_value) + " | " + str(value))


cpdef double[:] recalc_common_co_oldnew(int[:] new_starts, long[:] new_hash_store, int[:] old_starts, long[:] old_hash_store):
    global co_matrix
    cdef int i, j, n, m
    n = <int> len(new_starts)
    m = <int> len(old_starts)
    cdef double[:] result = init_double_array(n*m-1+1, 0.0)

    with nogil, parallel(num_threads=12):
        global co_matrix
        n = <int> len(new_starts)
        m = <int> len(old_starts)
        #for i in prange(0, n, nogil=False, schedule='guided', chunksize=25):
        for i in range(0, n):
            #cluster1 = clip_hashes(i, starts, hash_store)
            for j in range(m):
            #   cluster2 = clip_hashes(j, starts, hash_store)
                result[i*n+j] = CalcCommonCO(\
                    new_hash_store[ new_starts[i]:get_end_index(new_starts, i) ],\
                    old_hash_store[ old_starts[j]:get_end_index(old_starts, j) ],\
                    co_matrix)
            with gil:
                print("finished " + str(i))
    return result

cpdef double[:] recalc_common_co_newnew(int[:] starts, long[:] hash_store):
    global co_matrix
    cdef int i, j, n, size, index
    n = len(starts)
    size = floor((n*(n+1))/2)
    cdef double[:] result = init_double_array(size, 0.0)

    for i in prange(n, nogil=True):
        n = <int> len(starts)
        #cluster1 = clip_hashes(i, starts, hash_store)
        for j in range(i+1,  n):
        #   cluster2 = clip_hashes(j, starts, hash_store)
            index = map_index(i,j,n)
            result[index] = CalcCommonCO(hash_store[ starts[i]:get_end_index(starts, i) ],\
                hash_store[ starts[j]:get_end_index(starts, j)], co_matrix)
    
    return result

cdef int get_end_index(int[:] starts, int i) nogil:
    cdef int n = <int> len(starts)
    return starts[i+1] if i < (n-1) else n


cdef int map_index(int x, int y, int size) nogil:
    cdef int i, j, row
    i = x if x < y else y
    j = y if x < y else x

    row = i*size
    
    return row - floor( ( i*(i+1) ) / 2) + j

cdef double CalcCommonCO(long[:] cluster1, long[:] cluster2, HashTable co_matrix) nogil:
    cdef int i, j, l1, l2
    cdef double value = 0.0
    cdef int index = 0
    l1 = <int> len(cluster1)
    l2 = <int> len(cluster2)

    if l1 == 0 or l2 == 0: 
        return 0.0

    for i in range(l1):
        for j in range(l2):
            index = hash(cluster1[i], cluster2[j])
            value += co_matrix.get(index, 0.0)
            
    return value / (len(cluster1) + len(cluster2))
