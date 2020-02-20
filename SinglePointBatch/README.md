# Single point batch

To run a batch of single point calculation:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Prepare an appended geometry file
3. Prepare a job script for your queuing system
4. Entre a work directory
5. `python3 create.py`
6. `bash submit.sh`

Optionally, after not every job has finished, remove unnecessary files and directories for collecting by `bash clean.sh`

After a batch of Columbus7 single point calculation, collect data by `python3 collect.py`. Single point information including energy, energy gradient, interstate coupling, transition dipole will be appended

Optionally, visualize the potential energy surface by `python3 surface.py`

Tips:
1. create.py creates directories naming from 0 to the number of geometries - 1, which are single point job directories
2. submit.sh submits unsubmitted job and failed job

Examples to prepare appended geometry file: 
* `python3 intscan.py` generates a scan along an internal coordinate
* `python3 path.py` generates a linear synchronous transit path from geometry 1 to geometry 2
