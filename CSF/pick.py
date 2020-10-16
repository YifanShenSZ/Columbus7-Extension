import argparse
from pathlib import Path
from typing import List

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("CSF_txt", type=Path, help="CSF list from the end of drtls")
    args = parser.parse_args()
    return args

def select(CSF:List) -> bool:
    # # N-H1 sigma -> pi*
    # if CSF[0] + CSF[7] + CSF[13] < 2: return True
    # # N-H2 sigma -> pi*
    # if CSF[1] + CSF[7] + CSF[12] < 2: return True
    # # C-N sigma -> pi*
    # if CSF[2] + CSF[7] + CSF[11] < 2: return True
    # # N-H1 and N-H2 sigma -> pi*
    # if CSF[0] + CSF[1] + CSF[7] + CSF[12] + CSF[13] < 4: return True
    # # N-H1 and C-N sigma -> pi*
    # if CSF[0] + CSF[2] + CSF[7] + CSF[11] + CSF[13] < 4: return True
    # # N-H2 and C-N sigma -> pi*
    # if CSF[1] + CSF[2] + CSF[7] + CSF[11] + CSF[12] < 4: return True
    # # N-H1 and N-H2 and C-N sigma -> pi*
    # if CSF[0] + CSF[1] + CSF[2] + CSF[7] + CSF[11] + CSF[12] + CSF[13] < 6: return True

    # N-H sigma -> pi*
    # if CSF[0] + CSF[1] + CSF[6] + CSF[10] + CSF[11] < 4: return True

    # N-H sigma --> pi*
    if CSF[0] + CSF[1] + CSF[6] + CSF[10] + CSF[11] < 3: return True
    return False

if __name__ == "__main__":
    args = parse_args()
    with open("CSF.txt", 'r') as f: lines = f.readlines()
    CSFs = []
    for line in lines:
        line = line.split()
        for i in range(2, len(line)):
            line[i] = int(line[i])
            if 0 < line[i] < 3: line[i] = 1
            if line[i] == 3: line[i] = 2
        CSFs.append(line[2:])

    for i in range(len(CSFs)):
        CSF = CSFs[i]
        if select(CSF): print(i + 1)