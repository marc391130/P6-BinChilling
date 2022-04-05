from typing import Dict, List
from math import sqrt, pow

class DistanceCalculator:
    def calc_distance(self, element1: List[float], element2: List[float]):
        pass
    
    def __assert_same_length__(self, element1: List, element2: List):
        if len(element1) != len(element2):
            raise Exception(f"The two lists are not of same length, element1-len: {len(element1)} element2-len: {len(element2)}")


class EuclideanCalculator(DistanceCalculator):
    def calc_distance(self, element1: List[float], element2: List[float]):
        self.__assert_same_length__(element1, element2)

        result_length = 0
        for i in range(len(element1)):
            result_length += pow(element1[i] - element2[i], 2)

        return sqrt(result_length)

class CosignCalculator(DistanceCalculator):
    def calc_distance(self, element1: List[float], element2: List[float]):
        self.__assert_same_length__(element1, element2)

        counter = 0
        denominatorX = 0
        denominatorY = 0
        for i in range(len(element1)):
            counter += element1[i] * element2[i]
            denominatorX += pow(element1[i], 2)
            denominatorY += pow(element2[i], 2)

        return counter / (sqrt(denominatorX) * sqrt(denominatorY))
