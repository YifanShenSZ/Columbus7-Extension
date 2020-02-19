"""
Analyze vibration based on the internal coordinate Hessian
"""

''' Library '''
import argparse
from pathlib import Path
import shutil
import numpy
import FortranLibrary as FL
import basic

''' Global variable '''
intdim = 0 # Internal coordinate dimension

''' Routine '''
def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('IntCDefFormat', type=str, help='internal coordinate definition format: Columbus7 or default')
    parser.add_argument('IntCDef', type=Path, help='internal coordinate definition file')
    parser.add_argument('geom',    type=Path, help='Columbus7 geometry file')
    parser.add_argument('hessian', type=Path, help='Columbus7 hessian file')
    parser.add_argument('-o', '--output', type=str, default='geom.log', help='output file (default = geom.log)')
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
    # Read Hessian
    hessian = basic.read_hessian(args.hessian,intdim,intcdef)
    ''' Do the job '''
    # Get B^T matrix
    cartdim=3*NAtoms
    BT = numpy.empty((cartdim,intdim)); q = numpy.empty(intdim)
    FL.WilsonBMatrixAndInternalCoordinateq(r, BT, q, cartdim, intdim)
    # Calculate internal coordinate vibration
    freq = numpy.empty(intdim); LT = numpy.empty((intdim,intdim)); LinvT = numpy.empty((intdim,intdim))    
    FL.WilsonGFMethod(hessian, BT, mass, freq, LT, LinvT, intdim, NAtoms)
    # Transform to Cartesian coordinate
    cartmodeT = numpy.empty((intdim,cartdim))
    FL.InternalMode2CartesianMode(LT, BT, cartmodeT, intdim, cartdim)
    # Output for visualization
    r = r / FL.AInAU
    freq = freq / FL.cm_1InAu
    FL.Avogadro_Vibration(NAtoms, symbol, r, intdim, freq, cartmodeT, FileName=args.output)