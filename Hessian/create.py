'''
Create single point job directories, naming by Columbus7 convention:
CALC.c$[coordinate].d$[displacement]

Q: Why this script instead of colinp?
A: E.g. if we want some internal coordinate that is not supported by Columbus7
'''

from typing import List
import argparse
from pathlib import Path
import os

from numpy.lib.function_base import disp

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    # required arguments
    parser.add_argument("input", type=str, help="input template directory")
    parser.add_argument("IntCoordDef", type=str, help="internal coordinate definition file in default format")
    parser.add_argument("geom", type=str, help="reference xyz geometry")
    # optional arguments
    parser.add_argument("-d","--displacement", type=float, default=0.01, help="finite difference displacement (default = 0.01 a.u.)")
    parser.add_argument("-r","--rally_point", type=Path, default=Path("DISPLACEMENT"), help="location to gather the single points (default = DISPLACEMENT)")
    parser.add_argument("-s","--symmetric_modes", type=int, default=None, help="number of totally symmetric modes (default = assume C1)")
    # end
    args = parser.parse_args()
    return args

# read a vector from file (1 number/line)
def read_vector(file:Path) -> List:
    v = []
    with open(file,'r') as f:
        for line in f.readlines():
            v.append(float(line))
    return v
# write a vector to file (1 number/line)
def write_vector(v:List, file:Path) -> None:
    with open(file, 'w') as f:
        for element in v:
            print(element, file=f)

if __name__ == "__main__":
    ''' initialize '''
    args = parse_args()
    assert Path(args.input).exists()
    assert Path(args.IntCoordDef).exists()
    assert Path(args.geom).exists()
    os.system("~/Software/Mine/Tool-Collection/bin/cart2int.exe -f default -i " + args.IntCoordDef + " -x " + args.geom + " -o intgeom > cart2int.log")
    q_ref = read_vector("intgeom")
    intdim = q_ref.__len__()
    if args.symmetric_modes: symdim = args.symmetric_modes
    else:                    symdim = intdim
    if not args.rally_point.exists(): args.rally_point.mkdir()
    displfl = open(args.rally_point / "displfl", 'w')
    print(intdim, "/number of internal coordinates", sep=' ', file=displfl)
    print("no  /calculate dip.mom.derivatives", file=displfl)
    print("yes  /reference point calculation", file=displfl)
    ''' reference point '''
    point = args.rally_point / "REFPOINT"
    if not point.exists(): point.mkdir()
    point = str(point)
    command = "cp " + args.input + "/* " + point + "; " \
            + "cp " + args.IntCoordDef + " " + point + "; " \
            + "cp " + args.geom + " " + point + "; " \
            + "cd " + point + "; " \
            + "$COLUMBUS/xyz2col.x < " + args.geom
    os.system(command)
    ''' finite difference single points '''
    for i in range(symdim):
        q = q_ref.copy()
        # negative
        print(i + 1, -args.displacement, "doit", sep="   ", file=displfl)
        q[i] -= args.displacement
        point = args.rally_point / ("CALC.c" + str(i + 1) + ".d" + str(-args.displacement))
        if not point.exists(): point.mkdir()
        write_vector(q, point / "intgeom")
        point = str(point)
        command = "cp " + args.input + "/* " + point + "; " \
                + "cp " + args.IntCoordDef + " " + point + "; " \
                + "cp " + args.geom + " " + point + "; " \
                + "cd " + point + "; " \
                + "~/Software/Mine/Tool-Collection/bin/int2cart.exe -f default -i " + args.IntCoordDef + " -g intgeom -x " + args.geom + " -o cart.xyz > int2cart.log; " \
                + "$COLUMBUS/xyz2col.x < cart.xyz"
        os.system(command)
        # positive
        print(i + 1, args.displacement, "doit", sep="   ", file=displfl)
        q[i] += args.displacement + args.displacement
        point = args.rally_point / ("CALC.c" + str(i + 1) + ".d" + str(args.displacement))
        if not point.exists(): point.mkdir()
        write_vector(q, point / "intgeom")
        point = str(point)
        command = "cp " + args.input + "/* " + point + "; " \
                + "cp " + args.IntCoordDef + " " + point + "; " \
                + "cp " + args.geom + " " + point + "; " \
                + "cd " + point + "; " \
                + "~/Software/Mine/Tool-Collection/bin/int2cart.exe -f default -i " + args.IntCoordDef + " -g intgeom -x " + args.geom + " -o cart.xyz > int2cart.log; " \
                + "$COLUMBUS/xyz2col.x < cart.xyz"
        os.system(command)
    for i in range(symdim, intdim):
        q = q_ref.copy()
        # negative
        print(i + 1, -args.displacement, "doit", sep="   ", file=displfl)
        q[i] -= args.displacement
        point = args.rally_point / ("CALC.c" + str(i + 1) + ".d" + str(-args.displacement))
        if not point.exists(): point.mkdir()
        write_vector(q, point / "intgeom")
        point = str(point)
        command = "cp " + args.input + "/* " + point + "; " \
                + "cp " + args.IntCoordDef + " " + point + "; " \
                + "cp " + args.geom + " " + point + "; " \
                + "cd " + point + "; " \
                + "~/Software/Mine/Tool-Collection/bin/int2cart.exe -f default -i " + args.IntCoordDef + " -g intgeom -x " + args.geom + " -o cart.xyz > int2cart.log; " \
                + "$COLUMBUS/xyz2col.x < cart.xyz"
        os.system(command)
    displfl.close()
