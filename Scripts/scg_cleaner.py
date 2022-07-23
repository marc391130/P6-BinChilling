import sys
import os.path as path

def main(filepath):
    
    with open(filepath) as f:
        output = []
        for line in f:
            if line.startswith('edge'):
                left, right = line.split('\t')
                left = "edge_" + left.split('_')[1]
                output.append( left + '\t'+ right )
                
    outfile_name: str = path.basename(filepath)
    entities = outfile_name.split('.')
    outname = entities[0] + '.out' + (''.join(entities[1:]) if len(entities) > 1 else "" )
    with open(outname, 'w') as o:
        o.writelines(output)
                

if __name__ == "__main__":
    filepath = sys.argv[1]
    main(filepath)