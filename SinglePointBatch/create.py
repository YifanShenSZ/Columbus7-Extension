"""
Create single point job directories, 
naming from 0 to the number of geometries - 1 in current directory
"""

''' Library '''
import argparse
from pathlib import Path
import os

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('InputPath',type=Path,help='location of input template')
    parser.add_argument( 'GeomPath',type=Path,help='location of the appended geometry file')
    parser.add_argument('NAtoms',type=int,help='number of atoms in the molecule')
    args = parser.parse_args()
    return args

def main():
    # Initialize
    args = parse_args()
    with open(args.GeomPath,'r') as f: geoms = f.readlines()
    NGeom = int(len(geoms)/args.NAtoms)
    # Do the job
    line=0
    for i in range(NGeom):
        current=Path(str(i))
        if not current.exists(): current.mkdir()
        os.system('cp '+str(args.InputPath)+'/* '+str(current))
        with open(current/'geom','w') as f:
            for j in range(args.NAtoms):
                print(geoms[line],end='',file=f) # geoms[line] already ends with \n
                line=line+1

if __name__ == "__main__":
    main()