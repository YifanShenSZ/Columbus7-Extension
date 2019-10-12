# Columbus tricks

## Installation
### Standard step
1. tar -xzvf the installation package
2. export COLUMBUS=(current directory)/Columbus
3. Copy all GA libraries (in GA/lib) to $Columbus/source
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
2. For Rydberg basis, add Rydberg.bas to $COLUMBUS/source/iargos/basis

## Job control
1. Single point calculation: Never use NROOT to specified the number of states to compute, use FROOT instead (i.e. set both root and follow in colinp)
2. Geometry optimization: NROOT specifies which surface to optimize on. You can overwrite it with the transition moment input (i.e. set computing transition moment between n-th and n-th state)

## Computation cost
1. CI expansion should not exceed 100,000,000

## Historical issue
1. The number of each kind of orbitals appearing in cidrtin cannot exceed 256, including the frozen orbitals! (It's OK to code with unsigned short, but why you enforce it in every aspect??)

## Error message
1. When you see "drt size error" in ciudg execution, it is arised from slice > 6 or strtdrt + forbyt(bufsiz) - 1 > forbyt(drt%drtsiz)