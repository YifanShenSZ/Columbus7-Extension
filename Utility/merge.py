'''
Merge the loops along a path into a path, by extracting the first geometry in the loop
The merged data set is output to the parent directory of the loop directories
'''

import argparse
from pathlib import Path
import basic

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('LoopPath', type=Path, help='parent directory of the loop directories')
    parser.add_argument('StartDirectory', type=int, help='merge from directory')
    parser.add_argument('EndDirectory'  , type=int, help='StartDirectory to EndDirectory')
    parser.add_argument('NAtoms', type=int, help='number of atoms for a single geometry')
    parser.add_argument('-e', '--energy', action='store_true', help='only merge energy')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    # Command line input
    args = parse_args()
    ''' Do the job '''
    # geometry
    with open(args.LoopPath/"geom.data", 'a') as f:
        for i in range(args.StartDirectory, args.EndDirectory+1):
            with open(args.LoopPath/str(i)/"geom.data", 'r') as g:
                for j in range(args.NAtoms):
                    line = g.readline()
                    print(line, file=f, end='')
    # energy
    with open(args.LoopPath/"energy.data", 'a') as f:
        for i in range(args.StartDirectory, args.EndDirectory+1):
            with open(args.LoopPath/str(i)/"energy.data", 'r') as g:
                line = g.readline()
                print(line, file=f, end='')
    if not args.energy:
        NStates = len(line.split())
        # gradient
        for istate in range(NStates):
            path = 'cartgrad-'+str(istate+1)+'.data'
            with open(args.LoopPath/path,'a') as f:
                for i in range(args.StartDirectory, args.EndDirectory+1):
                    with open(args.LoopPath/str(i)/path, 'r') as g:
                        for j in range(args.NAtoms):
                            line = g.readline()
                            print(line, file=f, end='')
            for jstate in range(istate+1, NStates):
                path = 'cartgrad-'+str(istate+1)+'-'+str(jstate+1)+'.data'
                with open(args.LoopPath/path,'a') as f:
                    for i in range(args.StartDirectory, args.EndDirectory+1):
                        with open(args.LoopPath/str(i)/path, 'r') as g:
                            for j in range(args.NAtoms):
                                line = g.readline()
                                print(line, file=f, end='')
        # transition dipole
        for istate in range(NStates):
            for jstate in range(istate+1, NStates):
                path = 'transdip-'+str(istate+1)+'-'+str(jstate+1)+'.data'
                with open(args.LoopPath/path,'a') as f:
                    for i in range(args.StartDirectory, args.EndDirectory+1):
                        with open(args.LoopPath/str(i)/path, 'r') as g:
                            line = g.readline()
                            print(line, file=f, end='')