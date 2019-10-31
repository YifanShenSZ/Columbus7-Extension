# Columbus general experience

## Installation

* Steps
0. Install a Global Arrays
1. Enter a desired directory, tar -xzvf the installation package
2. export COLUMBUS=(current directory)/Columbus
3. cp PathToGlobalArrays/lib/*.a $Columbus
4. vi install.config. If the last line(s) = cpan || standard || grad || parallel, then delete these line(s). Repeat this operation each time before ./install.automatic
5. ./install.automatic cpan standard grad
6. ./install.automatic parallel

* Solutions to some issues
0. Global Arrays does not support mpi 3.0 standard, e.g. openmpi 4 fails but 3 works
1. Installation may fail with intel mpi. Openmpi works
2. Replace certain files with their counterparts in Modification

## Modification
The files in Modification directory is meant to replace their counterparts provided by Columbus in certain cases

* Installation
1. ./install.automatic may have trouble with late shell. Here is my modified version debugged on ubuntu 16 and 18
2. ./install.config is for intel mpi. Here is my version for open mpi
3. $COLUMBUS/machine.cfg/linux64.ifc.byterecl is a file to pass compiler flags. Replace it when some flag related problem occur, e.g. 'dynamic library has to be created with -fpic'

* Utilities
1. For minimum energy crossing point search, replace $COLUMBUS/polyhess.x
2. For special Rydberg basis, add Rydberg.bas to $COLUMBUS/source/iargos/basis. Columbus provides 2 common Rydberg basis in X.bas

## Weird stuff

* Basis and orbital
1. To use special basis, do not run prepinp
2. The number of each kind of orbitals appearing in cidrtin cannot exceed 256, including the frozen orbitals

* Geometry optimization
1. NROOT also specifies which surface to optimize on. You can overwrite it with the transition moment input transmomin (i.e. set computing transition moment between m-th and m-th state to optimize on m-th surface)