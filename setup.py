import subprocess
import os.path as path
import os
from time import sleep

SRC_FOLDER = path.join(os.getcwd(), 'src')
SLN_PATH = path.join(SRC_FOLDER, 'BinChillingTools.sln')
BC_FOLDER = path.join(os.getcwd(), 'BinChilling')
BUILD_FOLDER = path.join(SRC_FOLDER, 'Build')
EXECUTABLE_PATH = path.join(BC_FOLDER, 'BinChillingTools.exe')
PDB_PATH = path.join(BC_FOLDER, 'BinChillingTools.pdb')

def main():
    subprocess.run(['pip', 'install', '-r', './requirements.txt '])

    if path.exists(EXECUTABLE_PATH):
        print("Removing old EXE file: " + EXECUTABLE_PATH)
        os.remove(EXECUTABLE_PATH)
    if path.exists(PDB_PATH):
        print("Removing old PDB file: " + PDB_PATH)
        os.remove(PDB_PATH)
    if path.exists(BUILD_FOLDER) is False:
        print("Making build folder at: " + BUILD_FOLDER)
        subprocess.run(['mkdir', BUILD_FOLDER], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    args = ['dotnet', 'publish', SLN_PATH,
                    '-c', 'release', 
                    '-r', 'ubuntu.20.04-x64',
                    '-o', BUILD_FOLDER,
                    '-p:PublishSingleFile=true',
                    '-p:PublishTrimmed=true',
                    '--self-contained', 'true']
    print("Compiling using: " + ' '.join(args))
    subprocess.run(args)

    sleep(1)
    exe_output = path.join(BUILD_FOLDER, 'BinChillingTools')
    pdb_output = path.join(BUILD_FOLDER, 'BinChillingTools.pdb')
    os.rename(exe_output, EXECUTABLE_PATH)
    os.rename(pdb_output, PDB_PATH       )

    
    
if __name__ == '__main__':
    main()