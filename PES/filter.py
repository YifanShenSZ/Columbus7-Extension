'''
Filter the information in geom data to suit SurfGen format
'''

import argparse
from pathlib import Path
import basic

def parse_args() -> argparse.Namespace: # Command line input
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('geom_data', type=Path, help='geom data to filter')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    ''' Initialize '''
    # Command line input
    args = parse_args()
    ''' Do the job '''
    with open(args.geom_data, 'r') as f: lines = f.readlines()
    assert len(lines[0]) > 52, "This geom data file has already been filtered"
    with open(args.geom_data, 'w') as f:
        for line in lines:
            print(line[:2], line[10:52], sep='', file=f)