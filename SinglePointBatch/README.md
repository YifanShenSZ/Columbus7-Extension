# Single point batch

To run a batch of Columbus7 single point calculation:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Entre a work directory
3. cp Path_To_Columbus7-Extension/SinglePointBatch/{GenGeom.f90,makefile} .
4. Prepare intcfl (Columbus7 internal coordinate definition file) and geom (Columbus7 geometry file) in current directory. Modify GenGeom.f90 to generate a batch of desired geometries, make, ./GenGeom.exe (This step creates appended Columbus7 geometry file)
5. Entre a directory with enough free disk space
6. bash create.sh InputPath NAtoms GeomPath (This step creates directories naming from 1 to the number of geometries, which are Columbus7 single point job directories)
7. Prepare a job script for your queuing system, then bash submit.sh JobScriptPath (This step submits the batch jobs. If some jobs fail, bash again)

Optionally, after not every job has finished, you may bash clean.sh to remove unnecessary files and directories for collecting

All shell scripts support -h argument to show details

Dependency: my Fortran-Library, as written in makefile
