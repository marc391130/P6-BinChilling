# Auther: Narucs, Dak, Spand, & Trudsliv
# Created: 23-02-2022
# Intended purpose: To get data from fasta files
from ContigReader import ContigReader
import Constants as const
import numpy as np

cumreader = ContigReader()
# cumreader.save_as_numpy(const.FILEPATH, "test_data.npy")
# cumreader.load_numpy("test_data.npy")['edge_2'].pretty_print()
# print(cumreader.__get_abundance_dict__(const.ABUNDANCE_FILEPATH))

# r = FastaConverter()
# r.convert_data_to_new_file(const.FILEPATH, "temp_data.fasta")

d1 = cumreader.load_numpy('test_data_old.npy')
d2 = cumreader.load_numpy('test_data.npy')

if len(d1) != len(d2):
    raise Exception("d1 and d2 not same length")

with open("log.txt", 'x') as file:
    count = 0
    for key, value in d1.items():
        if d1[key].contig_length != d2[key].contig_length:
            file.write(f"Contig {key.rjust(10)} does not have the same ContigLen as its string. string length: {str(d1[key].contig_length).rjust(10)} edges_depth file ContigLen: {str(d2[key].contig_length).rjust(10)} differnce: {str((d1[key].contig_length - d2[key].contig_length)).rjust(10)}\n")
            print(f"{key} does not have the same length d1: {d1[key].contig_length:08d} d2: {d2[key].contig_length:08d} differnce: {(d1[key].contig_length - d2[key].contig_length):08d}")
            count += 1

    print(f"{count} failed out of {len(d1)}.")
    