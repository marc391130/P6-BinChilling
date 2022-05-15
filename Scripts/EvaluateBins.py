import sys
from typing import List, Dict, Tuple
from tqdm import tqdm
from CheckMFilter import BinDto, CheckMFilter
import argparse
import os

sys.path.insert(1, '../ACE')
from BinReaders import ContigReader, PartitionSetReader, SCGReader, ClusterReader
from ClusterDomain import Cluster
from Domain import ContigData
from BinEvaluator import BinEvaluator


def get_total_size(cluster) -> int:
    return sum([contig.contig_length for contig in cluster])

def main():
    parser = argparse.ArgumentParser( \
        prog='ACE', \
        description="""ACE BINNING ENSEMBLER""", \
        formatter_class=argparse.RawDescriptionHelpFormatter, \
        usage="%(prog)s WRITE THIS LATER PLZ", \
        add_help=True \
        )
    
    file_args = parser.add_argument_group(title='Contig input (required)', description=None)
    file_args.add_argument('--fasta', metavar='', required=True,\
        dest='fasta', help='path to fasta file of contigs')
    file_args.add_argument('--SCG', nargs='+', metavar='', required=True, \
        dest='SCG', help='Path to single copy genes file (required)')
    file_args.add_argument('--genefiles', metavar='', required=True, \
        dest='genefiles', help='Paths to gene files (required)')
    file_args.add_argument('--cache', metavar='', required=False, \
        dest='cache', help='Path to cache numpy file')
    file_args.add_argument('--cluster', metavar='', required=True, \
        dest='clusterpath', help='Path to cluster file (required)')
    file_args.add_argument('--depthfile', metavar='', required=True, \
        dest='depthfile', help='Path to depth file (required)')
    file_args.add_argument('--outputpath', metavar='', required=True, \
        dest='outputfile', help='Path to output file (required)')
    
    extra_args = parser.add_argument_group(title='Extra arguments', description=None)
    extra_args.add_argument('--minSize', metavar='', required=False,\
        dest='minSize', help='Minimum total size of clusters!')


    if len(sys.argv) <= 1:
        parser.print_help()
        return
    
    args = parser.parse_args()

    scg_files = [args.SCG] if isinstance(args.SCG, str) else args.SCG
    gene_db_files = [args.genefiles] if isinstance(args.genefiles, str) else args.genefiles
    paths = [args.fasta, args.depthfile, args.clusterpath] + scg_files + gene_db_files  
    print(paths)
    for path in paths:
        if os.path.isfile(path) is False:
            raise FileNotFoundError(path)

    minSize: int = 0
    if args.minSize is not None:
        minSize = int(args.minSize)

    scg_reader = SCGReader(scg_files, gene_db_files)
    reader = ContigReader(fasta_file=args.fasta, scg_reader=scg_reader, depth_file=args.depthfile, \
            numpy_file=args.cache)
    cluster_reader = ClusterReader(file_path=args.clusterpath, contig_reader=reader, numpy_file=args.cache)
    clusters = cluster_reader.clusters
    clusters = [cluster for cluster in clusters if get_total_size(cluster) >= minSize]
    scg_reader = reader.SCG_reader

    all_scgs = scg_reader.read_MS_scgs() # 
    evaluator = BinEvaluator(all_scgs)
    data = evaluator.evaluate_lst(clusters)

    result_dto_lst, total_com, total_con = [], 0, 0
    data_lst = list(data.items())
    for data_idx in range(len(data_lst)):
        cluster_idx, tuple = data_lst[data_idx]
        completeness, contamination, mp = tuple
        total_com, total_con = total_com + completeness,  total_con + contamination

        dto = BinDto(f"bin_{data_idx}", contamination, completeness)
        result_dto_lst.append(dto)
    
    checkMfilter = CheckMFilter(None, args.outputfile, lambda x: True)
    checkMfilter.write_output(result_dto_lst)
    print('\n')
    print(f'\nTotal completeness: {total_com} \nTotal contamination: {total_con}')


if __name__ == '__main__': 
    main()
    # if len(sys.argv) != 8:
    #     print("arguments need to be:\n", \
    #         "1: Fasta filepath\n", \
    #         "2: SCG_Filepath\n", \
    #         "3: gene db file\n",\
    #         "4: Numpy_filepath\n", \
    #         "5: Cluster_filepath\n", \
    #         "6: depth_filepath\n",\
    #         "7: Output path")
        
