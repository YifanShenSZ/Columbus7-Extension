"""
Analyze vibration based on the internal coordinate Hessian
"""

''' Library '''
import argparse
from pathlib import Path
import shutil
from typing import List
import numpy
import FortranLibrary as FL

''' Global variable '''
intdim = 0 # Internal coordinate dimension

''' Routine '''
def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('IntCDefFormat',type=str,help='internal coordinate definition format: Columbus7 or default')
    parser.add_argument('IntCDef',type=Path,help='internal coordinate definition file')
    parser.add_argument('geom',type=Path,help='Columbus7 geometry file')
    parser.add_argument('hessian',type=Path,help='Columbus7 hessian file')
    parser.add_argument('-o','--output',type=str,default='geom.log',help='output file (default = geom.log)')
    args = parser.parse_args()
    return args

def read_geom(GeomFile: Path) -> (int, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray):
    with open(GeomFile,'r') as f: lines=f.readlines()
    NAtoms = len(lines)
    symbol = numpy.empty(NAtoms,dtype=str)
    number = numpy.empty(NAtoms)
    r      = numpy.empty(int(3*NAtoms))
    mass   = numpy.empty(NAtoms)
    for i in range(NAtoms):
        temp = lines[i].split()
        symbol[i] = temp[0].strip()
        number[i] = float(temp[1].strip())
        r[3*i  ]  = float(temp[2].strip())
        r[3*i+1]  = float(temp[3].strip())
        r[3*i+2]  = float(temp[4].strip())
        mass[i]   = float(temp[5].strip())
    mass *= FL.AMUInAU
    return NAtoms, symbol, number, r, mass

def read_hessian(HessianFile: Path, intdim: int, intcdef: List) -> numpy.ndarray:
    hessian = numpy.empty((intdim,intdim))
    with open(HessianFile,'r') as f: lines=f.readlines()
    lineindex = 0
    for i in range(intdim):
        for j in range(0,intdim,8):
            temp = lines[lineindex].split(); lineindex += 1
            jstart=j; jstop=j+8
            if jstop>intdim: jstop=intdim
            index = 0
            for jj in range(jstart,jstop):
                hessian[i,jj] = float(temp[index].strip())
                index += 1
    # The internal coordinate and vibration routines of Columbus use weird unit:
    #     energy in 10^-18 J, length in A (to be continued)
    hessian /= 4.35974417 # 1 Hatree = 4.35974417 * 10^-18 J
    for i in range(intdim):
        if intcdef[i].motion[0].type == 'stretching':
            hessian[i,:] /= FL.AInAU
            hessian[:,i] /= FL.AInAU  
    return hessian

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
    NAtoms, symbol, number, r, mass = read_geom(args.geom)
    # Read Hessian
    hessian = read_hessian(args.hessian,intdim,intcdef)
    ''' Do the job '''
    # Get B^T matrix
    cartdim=3*NAtoms
    BT = numpy.empty((cartdim,intdim)); q = numpy.empty(intdim)
    FL.WilsonBMatrixAndInternalCoordinateq(r, BT, q, cartdim, intdim)
    # Calculate internal coordinate vibration
    freq = numpy.empty(intdim); LT = numpy.empty((intdim,intdim)); LinvT = numpy.empty((intdim,intdim))    
    FL.WilsonGFMethod(hessian, BT, mass, freq, LT, LinvT, intdim, NAtoms)
    
    print(freq / FL.cm_1InAu)
    
    # Transform to Cartesian coordinate
    cartmodeT = numpy.empty((intdim,cartdim))
    FL.InternalMode2CartesianMode(LT, BT, cartmodeT, intdim, cartdim)
    # Output for visualization
    r = r / FL.AInAU
    freq = freq / FL.cm_1InAu
    FL.Avogadro_Vibration(NAtoms, symbol, r, intdim, freq, cartmodeT, FileName=args.output)