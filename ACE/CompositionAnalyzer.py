import Constants as const
from Domain import Composition

class CompositionAnalyzer:
    def __init__(self) -> None:
        pass

    def analyze_composition(self, composition_dict : Composition, contig_string : str) -> None:
        com_len = len(contig_string) - (const.COMPOSITION_CONSTANT - 1)
        #com_normalize_val = 1 / com_len
        if com_len < const.COMPOSITION_CONSTANT:
            raise Exception("composition string is smaller than the COMPOSITION_CONSTANT: " + str(const.COMPOSITION_CONSTANT))

        for x in range(com_len):
            key = contig_string[x:x + const.COMPOSITION_CONSTANT]
            reversed_key = key[::-1]
            if key in composition_dict:
                composition_dict.AddOccurence(key)
            elif reversed_key in composition_dict:
                composition_dict.AddOccurence(reversed_key)