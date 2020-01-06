"""
Collect Columbus7 result from the single point batch
Default mode: collect MRCI energy, gradient, transition dipole
Alternative modes are available, see optional arguments
Please note that interstate coupling replaces nonadiabatic coupling
"""

''' Library '''
import argparse
from pathlib import Path
import numpy

''' Global variable '''
NData  = 0 # Number of data points
NAtoms = 0 # Number of atoms
intdim = 0 # Internal coordinate dimension

''' Routine '''
def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('BatchPath',type=Path,help='location of the single point batch directory')
    parser.add_argument('StartDirectory',type=int,help='collect from directory')
    parser.add_argument('EndDirectory',type=int,help='StartDirectory to EndDirectory')
    parser.add_argument('NState',type=int,help='collect from state 1 to NState')
    parser.add_argument('-s','--single', action='store_true',help='only collect NState-th state')
    parser.add_argument('-e','--energy' ,action='store_true',help='only collect energy')
    parser.add_argument('-i','--intgrad',action='store_true',help='additionally collect internal coordinate gradient')
    parser.add_argument('-m','--mcscf',  action='store_true',help='collect MCSCF result rather than MRCI')
    args = parser.parse_args()
    return args

def SingleState(args: argparse.Namespace):
    # Allocate memory
    energy=numpy.empty(NData)
    cartgrad=numpy.empty((NData,NAtoms,3))
    if args.intgrad: intgrad=numpy.empty((NData,intdim))
    if args.energy: # only collect energy
        # Read
        for i in range(args.StartDirectory,args.EndDirectory+1):
            with open(args.BatchPath/str(i)/'LISTINGS'/'ciudgsm.sp','r') as f: lines=f.readlines()
            for j in range(len(lines)):
                if 'mr-sdci  convergence criteria satisfied' in lines[j]: break
            if j == len(lines) - 1 :
                print('Warning: mr-sdci did not converge at job %d' % i)
            else:
                temp=lines[j+2+args.NState].split()
                energy[i]=float(temp[len(temp)-5].replace('D','e'))
        # Output
        with open('energy.all','w') as f:
            for i in range(NData): print(energy[i],file=f)
    else: # collect energy and gradient
        # Read
        for i in range(args.StartDirectory,args.EndDirectory+1):
            # energy
            with open(args.BatchPath/str(i)/'LISTINGS'/'ciudgsm.sp','r') as f: lines=f.readlines()
            for j in range(len(lines)):
                if 'mr-sdci  convergence criteria satisfied' in lines[j]: break
            if j == len(lines) - 1 :
                print('Warning: mr-sdci did not converge at job %d' % i)
            else:
                temp=lines[j+2+args.NState].split()
                energy[i]=float(temp[len(temp)-5].replace('D','e'))
                # gradient
                with open(args.BatchPath/str(i)/'GRADIENTS'/('cartgrd.drt1.state'+str(args.NState)+'.sp'),'r') as f: lines=f.readlines()
                for j in range(NAtoms):
                    temp=lines[j].split()
                    for k in range(3): cartgrad[i,j,k]=float(temp[k].replace('D','e'))
                if args.intgrad: # internal coordinate gradient
                    with open(args.BatchPath/str(i)/'GRADIENTS'/('intgrd.drt1.state'+str(args.NState)+'.sp'),'r') as f: lines=f.readlines()
                    for j in range(intdim): intgrad[i,j]=float(lines[j].replace('D','e'))
        # Output
        with open('energy.all','w') as f: # energy
            for i in range(NData): print(energy[i],file=f)
        with open('cartgrd.drt1.state'+str(args.NState)+'.all','w') as f: # gradient
            for i in range(NData):
                for j in range(NAtoms): print(cartgrad[i,j,0],cartgrad[i,j,1],cartgrad[i,j,2],sep=' ',file=f)
        if args.intgrad: # internal coordinate gradient
            with open('intgrd.drt1.state'+str(args.NState)+'.all','w') as f:
                for i in range(NData):
                    for j in range(intdim): print(intgrad[i,j],file=f)

def MultiState(args: argparse.Namespace):
    # Allocate memory
    energy=numpy.empty((NData,args.NState))
    cartgrad=numpy.empty((NData,args.NState,args.NState,NAtoms,3))
    if args.intgrad: intgrad=numpy.empty((NData,args.NState,args.NState,intdim))
    if args.energy: # only collect energy
        # Read
        for i in range(args.StartDirectory,args.EndDirectory+1):
            with open(args.BatchPath/str(i)/'LISTINGS'/'ciudgsm.sp','r') as f: lines=f.readlines()
            for j in range(len(lines)):
                if 'mr-sdci  convergence criteria satisfied' in lines[j]: break
            if j == len(lines) - 1 :
                print('Warning: mr-sdci did not converge at job %d' % i)
            else:
                for k in range(args.NState):
                    temp=lines[j+2+k+1].split()
                    energy[i,k]=float(temp[len(temp)-5].replace('D','e'))
        # Output
        with open('energy.all','w') as f:
            for i in range(NData):
                for j in range(args.NState-1): print(energy[i,j],end=' ',file=f)
                print(energy[i,args.NState-1],file=f)
    else: # collect energy and gradient
        # Read
        for i in range(args.StartDirectory,args.EndDirectory+1):
            # energy
            with open(args.BatchPath/str(i)/'LISTINGS'/'ciudgsm.sp','r') as f: lines=f.readlines()
            for j in range(len(lines)):
                if 'mr-sdci  convergence criteria satisfied' in lines[j]: break
            if j == len(lines) - 1 :
                print('Warning: mr-sdci did not converge at job %d' % i)
            else:
                for k in range(args.NState):
                    temp=lines[j+2+k+1].split()
                    energy[i,k]=float(temp[len(temp)-5].replace('D','e'))
                # gradient
                for istate in range(args.NState):
                    with open(args.BatchPath/str(i)/'GRADIENTS'/('cartgrd.drt1.state'+str(istate+1)+'.sp'),'r') as f: lines=f.readlines()
                    for j in range(NAtoms):
                        temp=lines[j].split()
                        for k in range(3): cartgrad[i,istate,istate,j,k]=float(temp[k].replace('D','e'))
                    for jstate in range(istate+1,args.NState):
                        with open(args.BatchPath/str(i)/'GRADIENTS'/('cartgrd.nad.drt1.state'+str(istate+1)+'.drt1.state'+str(jstate+1)+'.sp'),'r') as f: lines=f.readlines()
                        for j in range(NAtoms):
                            temp=lines[j].split()
                            for k in range(3): cartgrad[i,istate,jstate,j,k]=float(temp[k].replace('D','e'))
                if args.intgrad: # internal coordinate gradient
                    for istate in range(args.NState):
                        with open(args.BatchPath/str(i)/'GRADIENTS'/('intgrd.drt1.state'+str(istate+1)+'.sp'),'r') as f: lines=f.readlines()
                        for j in range(intdim): intgrad[i,istate,istate,j]=float(lines[j].replace('D','e'))
                        for jstate in range(istate+1,args.NState):
                            with open(args.BatchPath/str(i)/'GRADIENTS'/('intgrd.nad.drt1.state'+str(istate+1)+'.drt1.state'+str(jstate+1)+'.sp'),'r') as f: lines=f.readlines()
                            for j in range(intdim): intgrad[i,istate,jstate,j]=float(lines[j].replace('D','e'))
        # Output
        with open('energy.all','w') as f: # energy
            for i in range(NData):
                for j in range(args.NState-1): print(energy[i,j],end=' ',file=f)
                print(energy[i,args.NState-1],file=f)
        for istate in range(args.NState): # gradient
            with open('cartgrd.drt1.state'+str(istate+1)+'.all','w') as f:
                for i in range(NData):
                    for j in range(NAtoms): print(cartgrad[i,istate,istate,j,0],cartgrad[i,istate,istate,j,1],cartgrad[i,istate,istate,j,2],sep=' ',file=f)
            for jstate in range(istate+1,args.NState):
                with open('cartgrd.nad.drt1.state'+str(istate+1)+'.drt1.state'+str(jstate+1)+'.all','w') as f:
                    for i in range(NData):
                        cartgrad[i,istate,jstate,:,:]=cartgrad[i,istate,jstate,:,:]*(energy[i,jstate]-energy[i,istate])
                        for j in range(NAtoms): print(cartgrad[i,istate,jstate,j,0],cartgrad[i,istate,jstate,j,1],cartgrad[i,istate,jstate,j,2],sep=' ',file=f)
        if args.intgrad: # internal coordinate gradient
            for istate in range(args.NState):
                with open('intgrd.drt1.state'+str(istate+1)+'.all','w') as f:
                    for i in range(NData):
                        for j in range(intdim): print(intgrad[i,istate,istate,j],file=f)
                for jstate in range(istate+1,args.NState):
                    with open('intgrd.nad.drt1.state'+str(istate+1)+'.drt1.state'+str(jstate+1)+'.all','w') as f:
                        for i in range(NData):
                            intgrad[i,istate,jstate,:,:]=intgrad[i,istate,jstate,:,:]*(energy[i,jstate]-energy[i,istate])
                            for j in range(intdim): print(intgrad[i,istate,jstate,j],file=f)

def mcscf(args: argparse.Namespace):
    # Allocate memory
    energy=numpy.empty((NData,args.NState))
    # Read energy
    for i in range(args.StartDirectory,args.EndDirectory+1):
        with open(args.BatchPath/str(i)/'LISTINGS'/'mcscfsm.sp','r') as f: lines=f.readlines()
        for j in range(len(lines)):
            if 'Individual total energies for all states' in lines[j]: break
        if j == len(lines) - 1 :
            print('Warning: mcscf did not converge at job %d' % i)
        else:
            for k in range(args.NState):
                energy[i,k]=float(lines[j+k+1][42:61].strip())
    # Output energy
    with open('energy.all','w') as f:
        for i in range(NData):
            for j in range(args.NState-1): print(energy[i,j],end='\t',file=f)
            print(energy[i,args.NState-1],file=f)

if __name__ == "__main__":
    # Initialize
    args = parse_args()
    NData = args.EndDirectory - args.StartDirectory + 1
    with open(args.BatchPath/str(args.StartDirectory)/'geom','r') as f: # Get NAtoms
        NAtoms=len(f.readlines())
    if args.intgrad: # Get intdim
        with open(args.BatchPath/str(args.StartDirectory)/'GRADIENTS'/('intgrd.drt1.state'+str(args.NState)+'.sp'),'r') as f:
            intdim=len(f.readlines())
    # Do the job
    if args.mcscf:
        mcscf(args)
    elif args.single or (args.NState == 1):
        SingleState(args)
    else:
        MultiState(args)
