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
    parser.add_argument('IntCDefFormat', type=str, help='internal coordinate definition format: Columbus7 or default')
    parser.add_argument('IntCDef', type=Path, help='internal coordinate definition file')
    parser.add_argument('geom',    type=Path, help='Columbus7 geometry file')
    parser.add_argument('IntC2scan', type=int, help='internal coordinate to scan')
    parser.add_argument('-b','--bidirection', action='store_true', help='scan bidirectionsally (default = positive only)')
    parser.add_argument('-n','--NSteps', type=int, default=10  , help='number of scan steps (default = 10)')
    parser.add_argument('-l','--length', type=int, default=0.01, help='step length (default = 0.01)')
    parser.add_argument('-o','--output', type=Path, default='geom.data', help='output file (default = geom.data), will append if already exists')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    args = parse_args()
    # Define internal coordinate
    if args.IntCDefFormat == 'Columbus7': IntCDef = Path('intcfl') 
    else: IntCDef = Path('InternalCoordinateDefinition')
    if args.IntCDef != IntCDef:
        shutil.copy(args.IntCDef, IntCDef)
        intdim, intcdef = FL.FetchInternalCoordinateDefinition(args.IntCDefFormat)
        FL.DefineInternalCoordinate('Columbus7')
        IntCDef.unlink()
    else:
        intdim, intcdef = FL.FetchInternalCoordinateDefinition(args.IntCDefFormat)
        FL.DefineInternalCoordinate('Columbus7')
    # Read geometry
    NAtoms, symbol, number, r, mass = basic.read_geom(args.geom)
    cartdim = 3 * NAtoms
    q = numpy.empty(intdim)
    FL.InternalCoordinateq(r, q, cartdim, intdim)
    ''' Do the job '''
    q1 = numpy.empty(q.shape); r1 = numpy.empty(r.shape)
    if args.bidirection:
        for i in range(-args.NSteps, 0):
            q1[:] = q[:]; q1[args.IntC2scan-1] += i * args.length
            FL.CartesianCoordinater(q1, r1, intdim, cartdim, uniquify='assimilate', mass=mass, r0=r)
    for i in range(1, args.NSteps+1):
        q1[:] = q[:]; q1[args.IntC2scan-1] += i * args.length
        FL.CartesianCoordinater(q1, r1, intdim, cartdim, uniquify='assimilate', mass=mass, r0=r)
        basic.write_geom(args.output, NAtoms, symbol, number, r, mass)