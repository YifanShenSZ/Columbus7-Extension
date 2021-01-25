'''
Generate a linear synchronous transit path from geometry 1 to geometry 2
'''

import argparse
from pathlib import Path
import numpy
import FortranLibrary as FL
import basic

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('geom1', type=Path, help='Columbus7 geometry file 1')
    parser.add_argument('geom2', type=Path, help='Columbus7 geometry file 2')
    parser.add_argument('-n','--NSteps', type=int, default=10, help='number of linear synchronous transit steps (default = 10)')
    parser.add_argument('-o','--output', type=Path, default='geom.data', help='output file (default = geom.data), will append if already exists')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    # Command line input
    args = parse_args()
    # Define internal coordinate
    intdim, _ = FL.DefineInternalCoordinate(args.IntCoordDefForm, file=args.IntCoordDefFile)
    # Read geometry
    NAtoms, symbol, number, r1, mass = basic.read_geom(args.geom1)
    NAtoms, symbol, number, r2, mass = basic.read_geom(args.geom2)
    ''' Do the job '''
    dr = (r2 - r1) / (args.NSteps + 1)
    for i in range(1, args.NSteps + 1):
        r = r1 + i * dr
        basic.write_geom(args.output, NAtoms, symbol, number, r, mass)
