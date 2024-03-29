'''
Convert Cartesian coordinate gradient to my internal coordinate,
then calculate internal coordinate Hessian

This is a finite difference calculation from gradients,
the Hessian matrix elements will be symmetrized by Hij = (Hij + Hji) / 2
futher symmetrization will be performed given irreducible

Optionally, collect geometry, MRCI energy, gradient, transition dipole instead of Hessian

Output to DISPLACEMENT/../LISTINGS
'''

import argparse
from pathlib import Path
from typing import List, Tuple
import numpy
import os

args          = 0  # Command line input
listings      = 0  # Output path
intdim        = 0  # Internal coordinate dimension
displacements = [] # displacements[i][j] is the displacement along i-th internal coordinate

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("DISPLACEMENT_path", type=Path, help="location of single point job directories")
    parser.add_argument("NState", type=int, help="calculate for state 1 to NState")
    parser.add_argument("IntCoordDef", type=str, help="internal coordinate definition file in default format")
    parser.add_argument("geom", type=str, help="reference xyz geometry")
    parser.add_argument("-i","--irreducible", type=int, nargs='+', help="number of internal coordinates per irreducible")
    parser.add_argument("-c","--collect", action='store_true', help="collect geometry, MRCI energy, gradient, transition dipole instead of Hessian")
    args = parser.parse_args()
    return args

# read Columbus7 Cartesian gradient, return internal gradient
def cart2int(dir:Path, state:int) -> numpy.ndarray:
    # preprocess Columbus7 Cartesian gradient: nowadays floating point numbers use 'e' instead of 'D'
    cartgrad = dir/"GRADIENTS"/("cartgrd.drt1.state" + str(state+1) + ".sp")
    with open(cartgrad, 'r') as f: lines = f.readlines()
    cartgrad = dir/"grad.cart"
    with open(cartgrad, 'w') as f:
        for line in lines:
            print(line.replace('D', 'e'), end='', file=f)
    # call cart2int.exe
    command = "cd " + str(dir) + "; "
    command += "~/Software/Mine/Tool-Collection/bin/cart2int.exe -f default " \
             + "-i " + args.IntCoordDef + " " \
             + "-x " + args.geom + " " \
             + "-g grad.cart > cart2int.log; "
    command += "cd ../.."
    os.system(command)
    # read internal gradient
    intgrad = numpy.empty(intdim)
    with open(dir/"grad.int", 'r') as f: lines = f.readlines()
    for i in range(intdim): intgrad[i] = float(lines[i])
    return intgrad

# calculate finite difference Hessian
def Hessian() -> None:
    hessian = numpy.empty((args.NState, intdim, intdim))
    # if some coordinate has only 1 displacement,
    # then reference point is required, so read it
    for displacement in displacements:
        if len(displacement) == 1:
            intgrad_ref = numpy.empty((args.NState, intdim))
            dir = Path(args.DISPLACEMENT_path/"REFPOINT")
            for istate in range(args.NState):
                intgrad_ref[istate, :] = cart2int(dir, istate)
            break
    # read internal coordinate gradients and calculate finite difference
    for idim in range(intdim):
        # this assumes gradient to be antisymmetric
        if len(displacements[idim]) == 1:
            dir = Path(args.DISPLACEMENT_path/('CALC.c' + str(idim+1) + '.d' + str(displacements[idim][0])))
            for istate in range(args.NState):
                intgrad = cart2int(dir, istate)
                hessian[istate, idim, :] = (intgrad - intgrad_ref[istate, :]) / displacements[idim][0]
        # normal finite difference
        elif len(displacements[idim]) == 2:
            dir1 = Path(args.DISPLACEMENT_path/('CALC.c' + str(idim+1) + '.d' + str(displacements[idim][0])))
            dir2 = Path(args.DISPLACEMENT_path/('CALC.c' + str(idim+1) + '.d' + str(displacements[idim][1])))
            for istate in range(args.NState):
                intgrad1 = cart2int(dir1, istate)
                intgrad2 = cart2int(dir2, istate)
                hessian[istate, idim, :] = (intgrad1 - intgrad2) / (displacements[idim][0] - displacements[idim][1])
        else:
            print("There should be 1 or 2 displacements along internal coordinate", idim, sep=' ')
    # symmetrize
    for i in range(intdim): 
        for j in range(i+1, intdim):
            hessian[:, i, j] = (hessian[:, i, j] + hessian[:, j, i]) / 2.0
            hessian[:, j, i] =  hessian[:, i, j]
    if args.irreducible != None:
        start = 0
        for NCoord in args.irreducible:
            stop = start + NCoord
            hessian[:, start:stop, :start] = 0.0
            hessian[:, start:stop, stop: ] = 0.0
            start = stop
    # output
    for istate in range(args.NState):
        with open(listings/('hessian'+str(istate+1)),'w') as f:
            for i in range(intdim):
                for j in range(intdim):
                    print('%25.10f' % hessian[istate,i,j],end='',file=f)
                print(file=f)

# read directory, return (geometry, energy, gradient, dipole)
# the return data are original strings
def read_directory(direcotry: Path) -> Tuple[List, List, List, List]:
    ## geometry
    with open(direcotry/'geom', 'r') as f: geom = f.readlines()
    ## energy
    with open(direcotry/'LISTINGS'/'ciudgsm.sp','r') as f: lines=f.readlines()
    for title_line in range(len(lines)):
        if 'mr-sdci  convergence criteria satisfied' in lines[title_line]: break
    if title_line == len(lines) - 1 :
        print('Warning: mr-sdci did not converge at job ' + str(direcotry))
    else:
        energy = []
        for k in range(args.NState):
            temp = lines[title_line + 2 + k + 1].split()
            energy.append(temp[len(temp) - 5])
        ## gradient
        gradient = []
        for i in range(args.NState): gradient.append(0)
        for _ in range(len(gradient)):
            gradient[_] = []
            for i in range(args.NState): gradient[_].append(0)
        for istate in range(args.NState):
            with open(direcotry/'GRADIENTS'/('cartgrd.drt1.state' + str(istate+1) + '.sp'), 'r') as f:
                gradient[istate][istate] = f.readlines()
            for jstate in range(istate + 1, args.NState):
                with open(direcotry/'GRADIENTS'/('cartgrd.nad.drt1.state' + str(istate+1) + '.drt1.state' + str(jstate+1) + '.sp'), 'r') as f:
                    gradient[istate][jstate] = f.readlines()
        # transition dipole
        dipole = []
        for i in range(args.NState): dipole.append(0)
        for _ in range(len(dipole)):
            dipole[_] = []
            for i in range(args.NState): dipole[_].append(0)
        for istate in range(args.NState):
            for jstate in range(istate + 1,args.NState):
                with open(direcotry/'LISTINGS'/('trncils.FROMdrt1.state' + str(istate+1) + 'TOdrt1.state' + str(jstate+1)), 'r') as f: lines = f.readlines()
                for title_line in range(len(lines)):
                    if 'Transition moment components' in lines[title_line]: break
                temp = lines[title_line + 6].split()
                dipole[istate][jstate] = temp[2:5]
    return geom, energy, gradient, dipole

# collect geometry, MRCI energy, gradient, transition dipole
def collect() -> None:
    geoms     = []
    energies  = []
    gradients = []
    dipoles   = []
    # Read REFPOINT
    geom, energy, gradient, dipole = read_directory(Path(args.DISPLACEMENT_path/'REFPOINT'))
    geoms    .append(geom    )
    energies .append(energy  )
    gradients.append(gradient)
    dipoles  .append(dipole  )
    # Read displacement points
    for idim in range(intdim):
        for displacement in displacements[idim]:
            geom, energy, gradient, dipole = read_directory(
                Path(args.DISPLACEMENT_path/('CALC.c' + str(idim+1) + '.d' + str(displacement)))
                )
            geoms    .append(geom    )
            energies .append(energy  )
            gradients.append(gradient)
            dipoles  .append(dipole  )
    # Output
    with open(listings/'geom.data','w') as f: # geometry
        for geom in geoms:
            for atom in geom:
                print(atom[:2], atom[10:52], sep='', end='\n', file=f)
    with open(listings/'energy.data','w') as f: # energy
        for energy in energies:
            for state in energy: print(state, end='    ', file=f)
            print(file=f)
    for istate in range(args.NState): # gradient
        with open(listings/('cartgrad-' + str(istate+1) + '.data'), 'w') as f:
            for gradient in gradients:
                for atom in gradient[istate][istate]:
                    print(atom.replace('D', 'e'), end='', file=f)
        for jstate in range(istate + 1, args.NState):
            with open(listings/('cartgrad-' + str(istate+1) + '-' + str(jstate+1) + '.data'), 'w') as f:
                for gradient in gradients:
                    for atom in gradient[istate][jstate]:
                        print(atom.replace('D', 'e'), end='', file=f)
    for istate in range(args.NState): # transition dipole
        for jstate in range(istate+1,args.NState):
            with open(listings/('transdip-'+str(istate+1)+'-'+str(jstate+1)+'.data'),'w') as f:
                for dipole in dipoles:
                    for component in dipole[istate][jstate]: print(component, end='    ', file=f)
                    print(file=f)

if __name__ == "__main__":
    ''' Initialize '''
    args = parse_args() # command line input
    # get listings
    listings = args.DISPLACEMENT_path/'..'/'LISTINGS'
    if not listings.exists(): listings.mkdir()
    # get displacements
    with open(args.DISPLACEMENT_path/'displfl','r') as f:
        lines = f.readlines()
        strs = lines[0].split(); intdim = int(strs[0]) # get intdim
        displacements = []
        current_coord = -1
        for line in lines[3:]:
            strs = line.split()
            coord = int(strs[0]) - 1
            if coord != current_coord:
                current_coord = coord
                displacements.append([])
            displacements[current_coord].append(float(strs[1]))
    ''' Do the job '''
    if args.collect:
        collect()
    else:
        Hessian()
