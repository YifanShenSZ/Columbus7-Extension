"""
Generate a scan along an internal coordinate normal mode
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
    parser.add_argument('IntCoordDefForm', type=str, help='internal coordinate definition format: Columbus7 or default')
    parser.add_argument('IntCoordDef', type=Path, help='internal coordinate definition file')
    parser.add_argument('geom',    type=Path, help='Columbus7 geometry file')
    parser.add_argument('hessian', type=Path, help='Columbus7 hessian file')
    parser.add_argument('coord2scan', type=int, help='normal mode to scan')
    parser.add_argument('-c','--cartscan',    action='store_true', help='alternatively scan along Cartesian coordinate normal mode')
    parser.add_argument('-b','--bidirection', action='store_true', help='scan bidirectionsally (default = positive only)')
    parser.add_argument('-n','--NSteps', type=int, default=10, help='number of scan steps (default = 10)')
    parser.add_argument('-l','--length', type=int, default=1 , help='step length (default = 1)')
    parser.add_argument('-o','--output', type=Path, default='geom.data', help='output file (default = geom.data), will append if already exists')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    # Command line input
    args = parse_args()
    # Define internal coordinate
    intdim, intcoorddef = FL.FetchInternalCoordinateDefinition(args.IntCoordDefForm, file=args.IntCoordDef)
    FL.DefineInternalCoordinate(args.IntCoordDefForm, file=args.IntCoordDef)
    # Read geometry
    NAtoms, symbol, number, r, mass = basic.read_geom(args.geom)
    # Read Hessian
    hessian = basic.read_hessian(args.hessian,intdim,intcoorddef)
    # Get B^T matrix
    cartdim=3*NAtoms
    BT = numpy.empty((cartdim,intdim)); q = numpy.empty(intdim)
    FL.WilsonBMatrixAndInternalCoordinate(r, BT, q)
    # Get normal mode
    freq = numpy.empty(intdim)
    intmodeT = numpy.empty((intdim,intdim))
    LinvT = numpy.empty((intdim,intdim))
    cartmodeT = numpy.empty((intdim,cartdim))
    FL.WilsonGFMethod(hessian, BT, mass, freq, intmodeT, LinvT, cartmodeT)
    ''' Do the job '''
    q1 = numpy.empty(q.shape); r1 = numpy.empty(r.shape)
    if args.cartscan:
        if args.bidirection:
            for i in range(args.NSteps, 0, -1):
                r1[:] = r[:] - i * args.length * cartmodeT[args.coord2scan-1, :]
                basic.write_geom(args.output, NAtoms, symbol, number, r1, mass)
        for i in range(1, args.NSteps+1):
            r1[:] = r[:] + i * args.length * cartmodeT[args.coord2scan-1, :]
            basic.write_geom(args.output, NAtoms, symbol, number, r1, mass)
    else:
        if args.bidirection:
            rall = numpy.empty((args.NSteps,r.shape[0]))
            q1[:] = q[:] - args.length * intmodeT[args.coord2scan-1, :]
            FL.CartesianCoordinate(q1, rall[0,:], r0=r)
            for i in range(1, args.NSteps):
                q1[:] = q[:] - (i+1) * args.length * intmodeT[args.coord2scan-1, :]
                FL.CartesianCoordinate(q1, rall[i,:], r0=rall[i-1,:])
            for i in range(args.NSteps): basic.write_geom(args.output, NAtoms, symbol, number, rall[args.NSteps-1-i,:], mass)
        rsave = r.copy()
        for i in range(1, args.NSteps+1):
            q1[:] = q[:] + i * args.length * intmodeT[args.coord2scan-1, :]
            FL.CartesianCoordinate(q1, r1, r0=rsave)
            basic.write_geom(args.output, NAtoms, symbol, number, r1, mass)
            rsave[:] = r1[:]