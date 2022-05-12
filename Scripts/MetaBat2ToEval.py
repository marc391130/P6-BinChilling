import os
import os.path as path
import sys


def main():
    folder = sys.argv[1]
    outputfile = sys.argv[2]

    count = 0
    
    with open(outputfile, 'w') as o:
        for file in [x for x in os.listdir(folder) if x.startswith('.')]:
            count += 1
            with open(path.join(folder, file), 'r') as f:
                for line in f.readlines():
                    line = line.replace('\n', '')
                    o.write(f"{count}\t{line}\n")
    pass

if __name__ == "__main__":
    main()