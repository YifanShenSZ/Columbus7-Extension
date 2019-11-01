# Single point batch

To run a batch of Columbus7 single point calculation:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Prepare intcfl (Columbus7 internal coordinate definition file) and geom (Columbus7 geometry file) in your work directory. Modify GenGeom.f90 to generate a batch of desired geometries, make, ./GenGeom.exe (This step creates appended Columbus7 geometry file)
3. bash create InputPath NAtoms GeomPath (This step creates directories naming from 1 to the number of geometries, which are Columbus7 single point job directories)
4. Modify RunColumbus to fit your queue, then bash submit (This step submits the batch jobs. If some jobs fail, bash again)

Optionally, after not every job is finished, you may bash clean, so that unnecessary files and directories for collecting are deleted

Dependency: my Fortran-Library, as written in makefile