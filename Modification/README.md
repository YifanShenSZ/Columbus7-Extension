# Modification to Columbus7
The files here are meant to replace their Columbus7 counterparts in certain cases

Installation
* Default install.automatic may have trouble with late shell. Here is my modified version debugged on ubuntu 16 and 18
* Default install.config is for intel mpi. Here are my versions for openmpi with intel compiler

Utilities
* For special Rydberg basis, add Rydberg.bas to $COLUMBUS/source/iargos/basis after installation
* To tune AO-MO integral transformation cutoff threshold, replace $COLUMBUS/source/mcscf/mcscf2.f before installation (search 'YifanShenSZ' for my code)
* For minimum energy crossing point search, replace $COLUMBUS/polyhess.x after installation
