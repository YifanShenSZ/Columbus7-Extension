'''
Generate a scan along 1 or 2 internal direction(s)

If 1 direction, scan along positive or negative direction
If 2 directions, mesh along both negative and positive directions
'''

import argparse
from pathlib import Path
import numpy
import FortranLibrary as FL
import basic

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('IntCoordDefForm', type=str , help='internal coordinate definition format: Columbus7 or default')
    parser.add_argument('IntCoordDefFile', type=Path, help='internal coordinate definition file')
    parser.add_argument('geom', type=Path, help='Columbus7 geometry file')
    parser.add_argument('direction', type=Path, help='internal vector file')
    parser.add_argument('-n','--negative', action='store_true', help='scan along negative direction rather than positive')
    parser.add_argument('-N','--NSteps', type=int, default=10, help='number of scan steps (default = 10)')
    parser.add_argument('-l','--length', type=float, default=0.01, help='step length (default = 0.01)')
    parser.add_argument('-d2','--direction2', type=Path, help='internal vector file for direction 2')
    parser.add_argument('-n2','--NSteps2', type=int, default=10, help='number of scan steps for direction 2 (default = 10)')
    parser.add_argument('-l2','--length2', type=float, default=0.01, help='step length for direction 2 (default = 0.01)')
    parser.add_argument('-o','--output', type=Path, default='geom.data', help='output file (default = geom.data), will append if already exists')
    args = parser.parse_args()
    return args

def read_vector(VectorFile:Path) -> numpy.ndarray:
    with open(VectorFile,'r') as f: lines = f.readlines()
    n = len(lines)
    vector = numpy.empty(n)
    for i in range(n): vector[i] = float(lines[i].replace('D','e').strip())
    return vector

if __name__ == "__main__":
    ''' Initialize '''
    # Command line input
    args = parse_args()
    # Define internal coordinate
    intdim, _ = FL.DefineInternalCoordinate(args.IntCoordDefForm, file=args.IntCoordDefFile)
    # Read geometry
    NAtoms, symbol, number, r, mass = basic.read_geom(args.geom)
    cartdim = 3 * NAtoms
    q = numpy.empty(intdim)
    FL.InternalCoordinate(r, q)
    # Read direction
    d = read_vector(args.direction)
    d /= numpy.linalg.norm(d)
    if args.direction2 is not None:
        d2 = read_vector(args.direction2)
        d2 /= numpy.linalg.norm(d2)
    ''' Do the job '''
    q1 = numpy.empty(q.shape)
    r1 = numpy.empty(r.shape)
    if args.direction2 is None:
        rsave = r.copy()
        for i in range(1, args.NSteps+1):
            if args.negative:
                q1 = q - i * args.length * d
            else:
                q1 = q + i * args.length * d
            FL.CartesianCoordinate(q1, r1, r0=rsave)
            basic.write_geom(args.output, NAtoms, symbol, number, r1, mass)
            rsave[:] = r1[:]
    else:
        for i in range(-args.NSteps, args.NSteps+1):
            rsave = r.copy()
            for j in range(-args.NSteps2, args.NSteps2+1):
                print(i * args.length, j * args.length2, sep='\t')
                q1 = q + i * args.length * d + j * args.length2 * d2
                FL.CartesianCoordinate(q1, r1, r0=rsave)
                basic.write_geom(args.output, NAtoms, symbol, number, r1, mass)
                rsave[:] = r1[:]
                