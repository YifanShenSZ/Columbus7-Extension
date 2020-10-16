'''
Create single point job directories, 
naming from 0 to the number of geometries - 1 in current directory
'''

import argparse
from pathlib import Path
import os

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('input', type=Path, help='input template directory')
    parser.add_argument( 'geom', type=Path, help='appended geometry file')
    parser.add_argument('NAtoms', type=int, help='number of atoms in the molecule')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    # Initialize
    args = parse_args() # Command line input
    with open(args.geom,'r') as f: lines = f.readlines()
    NGeom = int(len(lines)/args.NAtoms)
    # Do the job
    line=0
    for i in range(NGeom):
        current=Path(str(i))
        if not current.exists(): current.mkdir()
        os.system('cp '+str(args.input)+'/* '+str(current))
        with open(current/'geom','w') as f:
            for j in range(args.NAtoms):
                print(lines[line],end='',file=f) # lines[line] already ends with \n
                line=line+1
