from typing import Generic, TypeVar, List
T = TypeVar('T')

class something(Generic[T]):
    def __init__(self, value: T):
        self.value = value

smt = something[str]("hello!")
smt2 = something[int](10)

if type(smt) == something:
    print(smt.value)
    
if type(smt2) == something:
    print(smt2.value)

if type(smt) == list:
    print('how')