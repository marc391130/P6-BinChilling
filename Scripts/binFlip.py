import sys


inputFile = sys.argv[1]
outputFile = sys.argv[2]

result = []
with open(outputFile, 'w') as o:
    with open(inputFile, 'r') as f:
        for line in f.readlines():
            temp = line.replace('\n', '').split('\t')
            o.write(f"{temp[1]}\t{temp[0]}\n")





