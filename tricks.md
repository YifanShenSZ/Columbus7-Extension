# Columbus tricks

## Installation
### Step
1. Enter a desired directory, tar -xzvf the installation package
2. export COLUMBUS=(current directory)/Columbus
3. Copy all GA libraries (in GA/lib) to $Columbus
4. vi install.config. If the last lines = cpan || standar || grad || parallel, then delete these lines. Repeat this operation each time before ./install.automatic
5. ./install.automatic cpan
6. ./install.automatic standard
7. ./install.automatic grad
8. ./install.automatic parallel

### Suggestion
1. Better use openmpi rather than intel mpi
2. The install.config packed in the installation package is for gnu compiler, refer to the version provided here in Modification for intel compiler

## Modification
1. For minimum energy crossing point search, replace the polyhess.x in $COLUMBUS
2. For special Rydberg basis, add Rydberg.bas to $COLUMBUS/source/iargos/basis. Columbus provides 2 common Rydberg basis in X.bas

## Weird bug
1. The number of each kind of orbitals appearing in cidrtin cannot exceed 256, including the frozen orbitals!
2. To use special basis, do not run prepinp
3. Single point calculation: Never use NROOT to specified the number of states to compute, use FROOT instead (i.e. set both root and follow in colinp)
4. Geometry optimization: NROOT specifies which surface to optimize on. You can overwrite it with the transition moment input (i.e. set computing transition moment between n-th and n-th state)

## Computation cost
1. CI expansion should not exceed 100,000,000