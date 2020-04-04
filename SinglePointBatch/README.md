# Single point batch

To run a batch of single point calculation:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Prepare an appended geometry file (e.g. see ../Example)
3. Prepare a job script for your queuing system
4. Entre a work directory
5. `python3 create.py`
6. `bash submit.sh`

Optionally, after not every job has finished, remove unnecessary files and directories for collecting by `bash clean.sh`

After a batch of Columbus7 single point calculation, collect data by `python3 collect.py`. Single point information including energy, energy gradient, interstate coupling, transition dipole will be appended

Tips:
1. create.py creates directories naming from 0 to the number of geometries - 1, which are single point job directories
2. submit.sh submits new job and failed job
