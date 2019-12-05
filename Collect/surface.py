"""
A quick potential energy surface plot
Input format = the output of collect.py (energy.all)
"""

''' Library '''
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Command line input
parser = argparse.ArgumentParser(__doc__)
parser.add_argument('DataPath',type=Path,help='location of the surface data')
args = parser.parse_args()

# Read energy
with open(args.DataPath,'r') as f: data=f.readlines()
n=len(data); NState=len(data[0].split())
x=np.empty(n); e=np.empty((NState,n))
for i in range(n):
    temp=data[i].split()
    x[i]=i
    for j in range(NState): e[j,i]=float(temp[j].strip())

# Plot
plt.figure()
for i in range(NState):
    plt.scatter(x,e[i,:])
plt.show()