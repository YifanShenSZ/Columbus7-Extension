# Columbus general experience

## Installation

* Step
1. Enter a desired directory, tar -xzvf the installation package
2. export COLUMBUS=(current directory)/Columbus
3. Copy all GA libraries (in GA/lib) to $Columbus
4. vi install.config. If the last lines = cpan || standard || grad || parallel, then delete these lines. Repeat this operation each time before ./install.automatic
5. ./install.automatic cpan
6. ./install.automatic standard
7. ./install.automatic grad
8. ./install.automatic parallel

* Suggestion
1. Better use openmpi rather than intel mpi
2. The install.config packed in the installation package is for gnu compiler, refer to the version provided here in Modification for intel compiler

## Modification
* For minimum energy crossing point search, replace the polyhess.x in $COLUMBUS
* For special Rydberg basis, add Rydberg.bas to $COLUMBUS/source/iargos/basis. Columbus provides 2 common Rydberg basis in X.bas

## Weird stuff

* Basis and orbital
1. To use special basis, do not run prepinp
2. The number of each kind of orbitals appearing in cidrtin cannot exceed 256, including the frozen orbitals

* Geometry optimization
1. NROOT also specifies which surface to optimize on. You can overwrite it with the transition moment input transmomin (i.e. set computing transition moment between m-th and m-th state to optimize on m-th surface)