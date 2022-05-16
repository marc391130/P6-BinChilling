from typing import Dict, List, Tuple
import argparse
import re

def parse_args() -> any:
    parser = argparse.ArgumentParser(
        prog='Metabinner SCG to SCG File!',
        description="""A Script for making SCG from metabinner SCG""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s WRITE THIS LATER PLZ",
        add_help=True
        )
    
    io_args = parser.add_argument_group(title='Required args', description=None)
    io_args.add_argument('--SCG', metavar='', required=True, type=str,\
        dest='SCG', help="The path to Metabinner SCG file!"
        )
    io_args.add_argument('--output', metavar='', required=True, type=str,\
        dest='output', help="The path to the output file!"
        )
    
    return parser.parse_args()



def make_output(args) -> None:
    def clean_name(name: str) -> str:
        return name.replace("'", '').replace("{", '').split('_')[0]

    contig_scg: Dict[str, set] = {}

    with open(args.SCG, 'r') as f:
        for line in f.readlines():
            # Here will be some string manipulation trickery!
            line = re.split("^[0-9]*\s", line)[1]

            segments_in_line = line.split('}, ')
            
            for segment in segments_in_line:
                segment_parts = re.split("(?<=[^0-9],\s)", segment)
                segment_start = segment_parts[0].split(': ')
                
                if segment_start[0] == '{}\n':
                    continue

                name, gene_set = clean_name(segment_start[0]), set([clean_name(segment_start[1])])

                for idx in range(1, len(segment_parts)):
                    gene = clean_name(segment_parts[idx].split(':')[0])
                    gene_set.add(gene)
                
                contig_scg[name] = contig_scg.get(name, set()).union(gene_set)
    
    with open(args.output, 'w') as f:
        for name, _set in contig_scg.items():
            for gene in _set:
                f.write(f"{name}\t{gene}\n")

if __name__ == '__main__':
    args = parse_args()
    make_output(args)
