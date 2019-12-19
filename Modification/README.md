# Modification to Columbus7
The files here are meant to replace their Columbus7 counterparts in certain cases

Installation
* Default install.automatic may have trouble with late shell. Here is my modified version debugged on ubuntu 16 and 18
* Default install.config is for intel mpi. Here are my versions for open mpi with intel compiler with avx and avx2 instruction sets
* Default $COLUMBUS/machine.cfg/linux64.ifc.byterecl is fine but may not be optimal. Here are my favourites

Utilities
* For special Rydberg basis, add Rydberg.bas to $COLUMBUS/source/iargos/basis after installation
* To tune ao-mo integral cutoff threshold, replace $COLUMBUS/source/mcscf/mcscf2.f before installation
* For minimum energy crossing point search, replace $COLUMBUS/polyhess.x after installation
