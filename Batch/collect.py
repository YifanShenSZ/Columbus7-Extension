'''
Collect Columbus7 result from the single point batch
Output to the parent directory of the single point batch directories
Default mode: collect geometry, MRCI energy, gradient, transition dipole
Alternative modes are available, see optional arguments
'''

import argparse
from pathlib import Path
from typing import List
import numpy

args   = 0  # Command line input
NData  = 0 # Number of data points
NAtoms = 0 # Number of atoms

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('BatchPath', type=Path, help='parent directory of the single point batch directories')
    parser.add_argument('StartDirectory', type=int, help='collect from directory')
    parser.add_argument('EndDirectory'  , type=int, help='StartDirectory to EndDirectory')
    parser.add_argument('NState', type=int, help='collect from state 1 to NState')
    parser.add_argument('-e','--energy', action='store_true', help='only collect energy')
    parser.add_argument('-m','--mcscf' , action='store_true', help='collect MCSCF result rather than MRCI')
    parser.add_argument('-s','--single', action='store_true',  help='NState-th state only')
    args = parser.parse_args()
    return args

# Read directory, return (geometry, energy, gradient, dipole)
# The return data are original strings
def read_directory(direcotry: Path) -> (List, List, List, List):
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

# Collect geometry, MRCI energy, gradient, transition dipole
def collect():
    # energy only
    if args.energy:
        energies = []
        # Read
        for i in range(args.StartDirectory, args.EndDirectory + 1):
            direcotry = args.BatchPath/str(i)
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
                energies.append(energy)
        # Output
        with open(args.BatchPath/'energy.data','w') as f:
            for energy in energies:
                for state in energy: print(state, end='    ', file=f)
                print(file=f)
    else: 
        geoms     = []
        energies  = []
        gradients = []
        dipoles   = []
        # Read
        for i in range(args.StartDirectory, args.EndDirectory + 1):
            geom, energy, gradient, dipole = read_directory(args.BatchPath/str(i))
            geoms    .append(geom    )
            energies .append(energy  )
            gradients.append(gradient)
            dipoles  .append(dipole  )
        # Output
        with open(args.BatchPath/'geom.data','w') as f: # geometry
            for geom in geoms:
                for atom in geom:
                    print(atom[:2], atom[10:52], sep='', end='\n', file=f)
        with open(args.BatchPath/'energy.data','w') as f: # energy
            for energy in energies:
                for state in energy: print(state, end='    ', file=f)
                print(file=f)
        for istate in range(args.NState): # gradient
            with open(args.BatchPath/('cartgrad-' + str(istate+1) + '.data'), 'w') as f:
                for gradient in gradients:
                    for atom in gradient[istate][istate]:
                        print(atom.replace('D', 'e'), end='', file=f)
            for jstate in range(istate + 1, args.NState):
                with open(args.BatchPath/('cartgrad-' + str(istate+1) + '-' + str(jstate+1) + '.data'), 'w') as f:
                    for gradient in gradients:
                        for atom in gradient[istate][jstate]:
                            print(atom.replace('D', 'e'), end='', file=f)
        for istate in range(args.NState): # transition dipole
            for jstate in range(istate+1,args.NState):
                with open(args.BatchPath/('transdip-'+str(istate+1)+'-'+str(jstate+1)+'.data'),'w') as f:
                    for dipole in dipoles:
                        for component in dipole[istate][jstate]: print(component, end='    ', file=f)
                        print(file=f)

# Read directory, return (geometry, energy, gradient)
# The return data are original strings
def read_directory_single(direcotry: Path) -> (List, str, str):
    # geometry
    with open(direcotry/'geom', 'r') as f: geom = f.readlines()
    # energy
    with open(direcotry/'LISTINGS'/'ciudgsm.sp','r') as f: lines=f.readlines()
    for title_line in range(len(lines)):
        if 'mr-sdci  convergence criteria satisfied' in lines[title_line]: break
    if title_line == len(lines) - 1 :
        print('Warning: mr-sdci did not converge at job ' + str(direcotry))
    else:
        temp = lines[title_line + 2 + args.NState].split()
        energy = temp[len(temp) - 5]
        # gradient
        with open(direcotry/'GRADIENTS'/('cartgrd.drt1.state' + str(args.NState) + '.sp'), 'r') as f:
            gradient = f.readlines()
    return geom, energy, gradient

# Collect geometry, MRCI energy, gradient
def collect_single() -> None:
    geoms     = []
    energies  = []
    gradients = []
    # Read
    for i in range(args.StartDirectory, args.EndDirectory + 1):
        geom, energy, gradient = read_directory_single(args.BatchPath/str(i))
        geoms    .append(geom    )
        energies .append(energy  )
        gradients.append(gradient)
    # Output
    # geometry
    with open(args.BatchPath/'geom.data','w') as f:
        for geom in geoms:
            for atom in geom:
                print(atom[:2], atom[10:52], sep='', end='\n', file=f)
    # energy
    with open(args.BatchPath/'energy.data','w') as f:
        for energy in energies:
            for state in energy: print(state, end='    ', file=f)
            print(file=f)
    # gradient
    with open(args.BatchPath/('cartgrad-' + str(args.NState) + '.data'), 'w') as f:
        for gradient in gradients:
            for atom in gradient:
                print(atom.replace('D', 'e'), end='', file=f)

# Currently, energy only
def mcscf(args: argparse.Namespace):
    energies = []
    # Read
    for i in range(args.StartDirectory, args.EndDirectory+1):
        with open(args.BatchPath/str(i)/'LISTINGS'/'mcscfsm.sp','r') as f: lines=f.readlines()
        for title_line in range(len(lines)):
            if 'Individual total energies for all states' in lines[title_line]: break
        if title_line == len(lines) - 1 :
            print('Warning: mcscf did not converge at job %d' % i)
        else:
            energy = []
            for k in range(args.NState):
                energy.append(lines[title_line + k + 1][42:61].strip())
            energies.append(energy)
    # Output
    with open(args.BatchPath/'energy.data','w') as f:
        for energy in energies:
            for state in energy: print(state, end='    ', file=f)
            print(file=f)

if __name__ == "__main__":
    # Initialize
    args = parse_args() # Command line input
    NData = args.EndDirectory - args.StartDirectory + 1
    with open(args.BatchPath/str(args.StartDirectory)/'geom','r') as f: # Get NAtoms
        NAtoms=len(f.readlines())
    # Do the job
    if args.mcscf:
        mcscf()
    elif args.single or (args.NState == 1):
        collect_single()
    else:
        collect()
