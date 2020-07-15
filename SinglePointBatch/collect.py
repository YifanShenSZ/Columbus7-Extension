"""
Collect Columbus7 result from the single point batch
Output to the parent directory of the single point batch directories
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

''' Routine '''
def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('BatchPath', type=Path, help='parent directory of the single point batch directories')
    parser.add_argument('StartDirectory', type=int, help='collect from directory')
    parser.add_argument('EndDirectory'  , type=int, help='StartDirectory to EndDirectory')
    parser.add_argument('NState', type=int, help='collect from state 1 to NState')
    parser.add_argument('-s', '--single', action='store_true', help='only collect NState-th state')
    parser.add_argument('-e', '--energy', action='store_true', help='only collect energy')
    parser.add_argument('-m', '--mcscf' , action='store_true', help='collect MCSCF result rather than MRCI')
    args = parser.parse_args()
    return args

def SingleState(args: argparse.Namespace):
    # Allocate memory
    energy=numpy.empty(NData)
    cartgrad=numpy.empty((NData,NAtoms,3))
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
                energy[i-args.StartDirectory]=float(temp[len(temp)-5].replace('D','e'))
        # Output
        with open(args.BatchPath/'energy.data','w') as f:
            for i in range(NData): print('%25.15E'%energy[i],file=f)
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
                energy[i-args.StartDirectory]=float(temp[len(temp)-5].replace('D','e'))
                # gradient
                with open(args.BatchPath/str(i)/'GRADIENTS'/('cartgrd.drt1.state'+str(args.NState)+'.sp'),'r') as f: lines=f.readlines()
                for j in range(NAtoms):
                    temp=lines[j].split()
                    for k in range(3): cartgrad[i-args.StartDirectory,j,k]=float(temp[k].replace('D','e'))
        # Output
        with open(args.BatchPath/'energy.data','w') as f: # energy
            for i in range(NData): print('%25.15E'%energy[i],file=f)
        with open(args.BatchPath/('cartgrad-'+str(args.NState)+'.data'),'w') as f: # gradient
            for i in range(NData):
                for j in range(NAtoms): print('%25.15E%25.15E%25.15E'%(cartgrad[i,j,0],cartgrad[i,j,1],cartgrad[i,j,2]),file=f)

def MultiState(args: argparse.Namespace):
    # Allocate memory
    energy  =numpy.empty((NData,args.NState))
    cartgrad=numpy.empty((NData,args.NState,args.NState,NAtoms,3))
    dipole  =numpy.empty((NData,args.NState,args.NState,3))
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
                    energy[i-args.StartDirectory,k]=float(temp[len(temp)-5].replace('D','e'))
        # Output
        with open(args.BatchPath/'energy.data','w') as f:
            for i in range(NData):
                for j in range(args.NState-1): print('%25.15E'%energy[i,j],end='',file=f)
                print('%25.15E'%energy[i,args.NState-1],file=f)
    else: # collect energy, gradient, transition dipole
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
                    energy[i-args.StartDirectory,k]=float(temp[len(temp)-5].replace('D','e'))
                # gradient
                for istate in range(args.NState):
                    with open(args.BatchPath/str(i)/'GRADIENTS'/('cartgrd.drt1.state'+str(istate+1)+'.sp'),'r') as f: lines=f.readlines()
                    for j in range(NAtoms):
                        temp=lines[j].split()
                        for k in range(3): cartgrad[i-args.StartDirectory,istate,istate,j,k]=float(temp[k].replace('D','e'))
                    for jstate in range(istate+1,args.NState):
                        with open(args.BatchPath/str(i)/'GRADIENTS'/('cartgrd.nad.drt1.state'+str(istate+1)+'.drt1.state'+str(jstate+1)+'.sp'),'r') as f: lines=f.readlines()
                        for j in range(NAtoms):
                            temp=lines[j].split()
                            for k in range(3): cartgrad[i-args.StartDirectory,istate,jstate,j,k]=float(temp[k].replace('D','e'))
                # transition dipole
                for istate in range(args.NState):
                    for jstate in range(istate+1,args.NState):
                        with open(args.BatchPath/str(i)/'LISTINGS'/('trncils.FROMdrt1.state'+str(istate+1)+'TOdrt1.state'+str(jstate+1)),'r') as f: lines=f.readlines()
                        for j in range(len(lines)):
                            if 'Transition moment components' in lines[j]: break
                        temp=lines[j+6].split()
                        dipole[i-args.StartDirectory,istate,jstate,0]=float(temp[2].replace('D','e'))
                        dipole[i-args.StartDirectory,istate,jstate,1]=float(temp[3].replace('D','e'))
                        dipole[i-args.StartDirectory,istate,jstate,2]=float(temp[4].replace('D','e'))
        # Output
        with open(args.BatchPath/'energy.data','w') as f: # energy
            for i in range(NData):
                for j in range(args.NState-1): print('%25.15E'%energy[i,j],end='',file=f)
                print('%25.15E'%energy[i,args.NState-1],file=f)
        for istate in range(args.NState): # gradient
            with open(args.BatchPath/('cartgrad-'+str(istate+1)+'.data'),'w') as f:
                for i in range(NData):
                    for j in range(NAtoms): print('%25.15E%25.15E%25.15E'%(cartgrad[i,istate,istate,j,0],cartgrad[i,istate,istate,j,1],cartgrad[i,istate,istate,j,2]),file=f)
            for jstate in range(istate+1,args.NState):
                with open(args.BatchPath/('cartgrad-'+str(istate+1)+'-'+str(jstate+1)+'.data'),'w') as f:
                    for i in range(NData):
                        cartgrad[i,istate,jstate,:,:]=cartgrad[i,istate,jstate,:,:]*(energy[i,jstate]-energy[i,istate])
                        for j in range(NAtoms): print('%25.15E%25.15E%25.15E'%(cartgrad[i,istate,jstate,j,0],cartgrad[i,istate,jstate,j,1],cartgrad[i,istate,jstate,j,2]),file=f)
        for istate in range(args.NState): # transition dipole
            for jstate in range(istate+1,args.NState):
                with open(args.BatchPath/('transdip-'+str(istate+1)+'-'+str(jstate+1)+'.data'),'w') as f:
                    for i in range(NData):
                        print('%25.15E%25.15E%25.15E'%(dipole[i,istate,jstate,0],dipole[i,istate,jstate,1],dipole[i,istate,jstate,2]),file=f)

def mcscf(args: argparse.Namespace): # Currently, energy only
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
                energy[i-args.StartDirectory,k]=float(lines[j+k+1][42:61].strip())
    # Output energy
    with open(args.BatchPath/'energy.data','w') as f:
        for i in range(NData):
            for j in range(args.NState-1): print('%25.15E'%energy[i,j],end='',file=f)
            print('%25.15E'%energy[i,args.NState-1],file=f)

if __name__ == "__main__":
    # Initialize
    args = parse_args() # Command line input
    NData = args.EndDirectory - args.StartDirectory + 1
    with open(args.BatchPath/str(args.StartDirectory)/'geom','r') as f: # Get NAtoms
        NAtoms=len(f.readlines())
    # Do the job
    if args.mcscf:
        mcscf(args)
    elif args.single or (args.NState == 1):
        SingleState(args)
    else:
        MultiState(args)
