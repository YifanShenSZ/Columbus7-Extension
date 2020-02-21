"""
Analyze vibration based on the internal coordinate Hessian
Return a file with geometry + frequencies + normal modes visualizable in Avogadro
"""

''' Library '''
import argparse
from pathlib import Path
import numpy
import FortranLibrary as FL
import basic

''' Global variable '''
intdim = 0 # Internal coordinate dimension

''' Routine '''
def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('IntCoordDefForm', type=str , help='internal coordinate definition format: Columbus7 or default')
    parser.add_argument('IntCoordDefFile', type=Path, help='internal coordinate definition file')
    parser.add_argument('geom',    type=Path, help='Columbus7 geometry file')
    parser.add_argument('hessian', type=Path, help='Columbus7 hessian file')
    parser.add_argument('-o', '--output', type=Path, default='geom.log', help='output file (default = geom.log)')
    parser.add_argument('-i', '--intcoord', action='store_true', help='additionally output internal coordinate normal modes')
    parser.add_argument('-io', '--intcoordoutput', type=Path, default='geom.int', help='internal coordinate normal modes output file (default = geom.int)')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    # Command line input
    args = parse_args()
    # Define internal coordinate
    intdim, intcdef = FL.FetchInternalCoordinateDefinition(args.IntCoordDefForm, file=args.IntCoordDefFile)
    FL.DefineInternalCoordinate(args.IntCoordDefForm, file=args.IntCoordDefFile)
    # Read geometry
    NAtoms, symbol, number, r, mass = basic.read_geom(args.geom)
    # Read Hessian
    hessian = basic.read_hessian(args.hessian,intdim,intcdef)
    ''' Do the job '''
    # Get B^T matrix
    cartdim=3*NAtoms
    BT = numpy.empty((cartdim,intdim)); q = numpy.empty(intdim)
    FL.WilsonBMatrixAndInternalCoordinateq(r, BT, q, cartdim, intdim)
    # Calculate vibration
    freq = numpy.empty(intdim)
    intmodeT = numpy.empty((intdim,intdim))
    LinvT = numpy.empty((intdim,intdim))
    cartmodeT = numpy.empty((intdim,cartdim))
    FL.WilsonGFMethod(hessian, BT, mass, freq, intmodeT, LinvT, cartmodeT, intdim, NAtoms)
    ''' Output '''
    # Convert to human unit
    r = r / FL.AInAU
    freq = freq / FL.cm_1InAu
    # Wilson GF method normalizes Cartesian coordinate normal mode by Hessian metric
    # However, this may not be an appropriate magnitude to visualize
    # Here we use infinity-norm to normalize Cartesian coordinate normal mode
    # Actually, normalize to 9.99 since the visualization file format is %5.2f
    for i in range(intdim): cartmodeT[i,:] *= 9.99 / numpy.amax(numpy.abs(cartmodeT[i,:]))
    FL.Avogadro_Vibration(NAtoms, symbol, r, intdim, freq, cartmodeT, file=args.output)
    if args.intcoord:
        # Wilson GF method normalizes internal coordinate normal mode by Hessian metric
        # However, this is inconvienient to tell the contribution of each internal coordinate
        # Here we use 1-norm to normalize internal coordinate normal mode
        temp, length = FL.dScientificNotation(float(intdim)); length += 1
        for i in range(intdim):
            temp = numpy.abs(intmodeT[i,:])
            intmodeT[i,:] /= numpy.sum(temp)
            # Output the dominating internal coordinate of each normal mode
            indextemp = numpy.argmax(temp)
            print('Normal mode %s is dominated by internal coordinate %s with %4.2f'%\
                (str(i+1).rjust(length), str(indextemp+1).rjust(length), numpy.abs(intmodeT[i,indextemp])))
        # Output a detailed contribution analysis to file
        with open(args.intcoordoutput,'w') as f:
            print('coord\\mode',end='\t',file=f)
            for i in range(1,intdim): print('%10i'%i,end='\t',file=f)
            print('%10i'%intdim,file=f)
            for i in range(intdim):
                print('%10i'%(i+1),end='\t',file=f)
                for j in range(intdim-1):
                    print('%10.7f'%intmodeT[j,i],end='\t',file=f)
                print('%10.7f'%intmodeT[intdim-1,i],file=f)