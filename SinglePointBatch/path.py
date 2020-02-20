"""
Generate a linear synchronous transit path from geometry 1 to geometry 2
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
    parser.add_argument('geom1',   type=Path, help='Columbus7 geometry file 1')
    parser.add_argument('geom2',   type=Path, help='Columbus7 geometry file 2')
    parser.add_argument('-n','--NSteps', type=int, default=10, help='number of linear synchronous transit steps (default = 10)')
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
        intdim = FL.DefineInternalCoordinate(args.IntCDefFormat)
        IntCDef.unlink()
    else:
        intdim = FL.DefineInternalCoordinate(args.IntCDefFormat)
    # Read geometry
    NAtoms, symbol, number, r1, mass = basic.read_geom(args.geom1)
    NAtoms, symbol, number, r2, mass = basic.read_geom(args.geom2)
    cartdim = 3 * NAtoms
    q1 = numpy.empty(intdim); FL.InternalCoordinateq(r1, q1, cartdim, intdim)
    q2 = numpy.empty(intdim); FL.InternalCoordinateq(r2, q2, cartdim, intdim)
    ''' Do the job '''
    dq = numpy.empty(intdim); dq[:] = (q2[:]-q1[:]) / (args.NSteps+1)
    q = numpy.empty(intdim)
    r = numpy.empty(cartdim); rsave = r1.copy()
    for i in range(1, args.NSteps+1):
        q[:] = q1[:] + float(i) * dq[:]
        FL.CartesianCoordinater(q, r, intdim, cartdim, uniquify='assimilate', mass=mass, r0=rsave)
        basic.write_geom(args.output, NAtoms, symbol, number, r, mass)
        rsave[:] = r[:]