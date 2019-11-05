# Modification to Columbus7
The files here are meant to replace their Columbus7 counterparts in certain cases

* Installation
1. install.automatic may have trouble with late shell. Here is my modified version debugged on ubuntu 16 and 18
2. install.config is for intel mpi. Here is my version for open mpi
3. $COLUMBUS/machine.cfg/linux64.ifc.byterecl is a file to pass compiler flags. Replace it when some flag related problem occur

* Utilities
1. For minimum energy crossing point search, replace $COLUMBUS/polyhess.x
2. For special Rydberg basis, add Rydberg.bas to $COLUMBUS/source/iargos/basis
