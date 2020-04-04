"""
A quick potential energy surface plot
Input format = the output of collect.py (energy.all)
"""

''' Library '''
import argparse
from pathlib import Path
import numpy
import matplotlib.pyplot as plt

''' Routine '''
def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('data', type=Path, help='energy data file')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args() # Command line input
    # Read energy
    with open(args.data,'r') as f: lines=f.readlines()
    n=len(lines); NState=len(lines[0].split())
    x=numpy.empty(n); e=numpy.empty((NState,n))
    for i in range(n):
        temp=lines[i].split()
        x[i]=i
        for j in range(NState): e[j,i]=float(temp[j].strip())
    # Plot
    plt.figure()
    for i in range(NState):
        plt.scatter(x,e[i,:])
    plt.show()
