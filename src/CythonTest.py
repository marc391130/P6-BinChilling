from sys import getsizeof
import BinChillingTools as BCT
import numpy as np
import multiprocessing as mp
import ctypes as c

SIZE  = 1_000_000_000


for x, i in ( (x, i) for x, i in enumerate(range(100)) if x % 2 == 0):
    print( f"{x}, {i}" )

# print(BCT.calc(5, 7))

# dic = {(i, i+1): i+2 for i in range(0, SIZE, 3)}


# start = np.array([x for x in range(0,SIZE,2)], dtype=np.intc)
# arr = np.array([x for x in range(SIZE)], dtype=np.int_)

# print("building...")

# BCT.init_coassosiation(dic, SIZE)

# print("calcing 1...")

# print(BCT.recalc_common_co_newnew(start, arr))

# print("calcing 2...")

# print(BCT.recalc_common_co_oldnew(start, arr, start, arr))
# print("{0:b}".format(-10))
# print("{0:b}".format(c.c_long(-10).value))
# print("{0:b}".format(c.c_ulong(-10).value))

# print(c.c_ulong(-10).value)