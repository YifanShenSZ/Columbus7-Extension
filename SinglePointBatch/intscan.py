"""
Generate a scan along an internal coordinate
"""

''' Library '''
import argparse
from pathlib import Path
import shutil
import numpy
import FortranLibrary as FL
import basic

''' Routine '''
def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('IntCoordDefForm', type=str , help='internal coordinate definition format: Columbus7 or default')
    parser.add_argument('IntCoordDefFile', type=Path, help='internal coordinate definition file')
    parser.add_argument('geom', type=Path, help='Columbus7 geometry file')
    parser.add_argument('coord2scan', type=int, help='internal coordinate to scan')
    parser.add_argument('-b','--bidirection', action='store_true', help='scan bidirectionsally (default = positive only)')
    parser.add_argument('-n','--NSteps', type=int, default=10  , help='number of scan steps (default = 10)')
    parser.add_argument('-l','--length', type=int, default=0.01, help='step length (default = 0.01)')
    parser.add_argument('-o','--output', type=Path, default='geom.data', help='output file (default = geom.data), will append if already exists')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    # Command line input
    args = parse_args()
    # Define internal coordinate
    intdim = FL.DefineInternalCoordinate(args.IntCoordDefForm, file=args.IntCoordDefFile)
    # Read geometry
    NAtoms, symbol, number, r, mass = basic.read_geom(args.geom)
    cartdim = 3 * NAtoms
    q = numpy.empty(intdim)
    FL.InternalCoordinateq(r, q, cartdim, intdim)
    ''' Do the job '''
    q1 = numpy.empty(q.shape); r1 = numpy.empty(r.shape)
    if args.bidirection:
        rall = numpy.empty((args.NSteps,r.shape[0]))
        q1[:] = q[:]; q1[args.coord2scan-1] -= args.length
        FL.CartesianCoordinater(q1, rall[0,:], intdim, cartdim, uniquify='assimilate', mass=mass, r0=r)
        for i in range(1, args.NSteps):
            q1[:] = q[:]; q1[args.coord2scan-1] -= (i+1) * args.length
            FL.CartesianCoordinater(q1, rall[i,:], intdim, cartdim, uniquify='assimilate', mass=mass, r0=rall[i-1,:])
        for i in range(args.NSteps): basic.write_geom(args.output, NAtoms, symbol, number, rall[args.NSteps-1-i,:], mass)
    rsave = r.copy()
    for i in range(1, args.NSteps+1):
        q1[:] = q[:]; q1[args.coord2scan-1] += i * args.length
        FL.CartesianCoordinater(q1, r1, intdim, cartdim, uniquify='assimilate', mass=mass, r0=rsave)
        basic.write_geom(args.output, NAtoms, symbol, number, r1, mass)
        rsave[:] = r1[:]