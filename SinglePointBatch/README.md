# Single point batch

To run a batch of single point calculation:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Prepare an appended geometry file
3. Prepare a job script for your queuing system
4. Entre a work directory
5. python create.py InputPath NAtoms GeomPath
6. bash submit.sh JobScriptPath

Optionally, after not every job has finished, you may remove unnecessary files and directories for collecting by:

    bash clean.sh

Tips:
1. create.py creates directories naming from 0 to the number of geometries - 1, which are single point job directories
2. submit.sh submits unsubmitted job and failed job

All scripts support -h argument to show details

Two examples to prepare appended geometry file: 
* scan.f90 shows how to generate a scan around a certain geometry
* path.f90 shows how to generate a linear path between two geometries
