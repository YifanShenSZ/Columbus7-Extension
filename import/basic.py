"""
General basic routines for Columbus7-Extension
"""

''' Library '''
from pathlib import Path
from typing import List
import numpy
import FortranLibrary as FL

''' Routine '''
# Read Columbus7 geometry file, return:
#     NAtoms (number of atoms)
#     symbol (element symbol of each atom)
#     number (element number of each atom)
#     r      (Cartesian coordinate in atomic unit)
#     mass   (mass of each atom in atomic unit)
def read_geom(GeomFile:Path) -> (int, numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray):
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

# Inverse to read_geom
def write_geom(GeomFile:Path, NAtoms:int, symbol:List, number:List, r:List, mass:List) -> None:
    with open(GeomFile,'a') as f:
        for i in range(NAtoms):
            print((' %-2s  %5.1f%14.8f%14.8f%14.8f%14.8f')%\
                (symbol[i],number[i],r[3*i],r[3*i+1],r[3*i+2],mass[i]/FL.AMUInAU),file=f)

# Read Columbus7 Hessian file, return Hessian in atomic unit
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
