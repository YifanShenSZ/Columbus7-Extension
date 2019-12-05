# Single point batch

To run a batch of Columbus7 single point calculation:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Prepare an appended Columbus7 geometry file
3. Prepare a job script for your queuing system
4. Entre a work directory
5. bash create.sh InputPath NAtoms GeomPath
6. bash submit.sh JobScriptPath

Optionally, after not every job has finished, you may remove unnecessary files and directories for collecting by:

    bash clean.sh

Tips:
1. create.sh creates directories naming from 1 to the number of geometries, which are Columbus7 single point job directories
2. submit.sh submits unsubmitted job and failed job

All shell scripts support -h argument to show details

Two examples to prepare appended Columbus7 geometry file: 
* loop.f90 shows how to generate a loop around a certain geometry
* path.f90 shows how to generate a linear path between two geometries
