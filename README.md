# Columbus7-Extension
Useful extension and modification for Columbus7

Featured utilities:
1. Run a batch of Columbus7 single point calculation (SinglePointBatch)
2. Collect the batch (Collect)
3. Calculate Hessian and vibration (Hessian), and provide scan along each normal mode

See each subdirectory for details

## Complement to Columbus7 documentation
Columbus7 originates in 1980 and preserves many historical codes. Besides, the user community of Columbus7 is always small. As a result, there are compatibility issues and undiscovered bugs, which are out of the coverage of the official documentation

### Installation
Steps
1. Install a Global Arrays
2. Enter a desired directory, tar -xzvf the installation package
3. export COLUMBUS=(current directory)/Columbus
4. cp Path_To_Global_Arrays/lib/*.a $Columbus
5. vi install.config. If the last line(s) = cpan || standard || grad || parallel, then delete these line(s). Repeat this operation each time before ./install.automatic
6. ./install.automatic cpan standard grad
7. ./install.automatic parallel

Solutions to some issues
* Global Arrays (subsequently Columbus7) does not support mpi 3.0 standard, e.g. openmpi 4 fails but 3 works
* Installation may fail with intel mpi. For now openmpi always works
* Replace certain files with their counterparts in Modification

### Basis set
To use special basis, do not run prepinp

The basis file ($COLUMBUS/source/iargos/basis/*.bas) format is:
1. Basis name
2. Number of basis types
3. Number of alphas, l + 1, number of contracted basis
4. Alpha value, contraction coefficients

The basis input daltaoin format is:
1. (s for spherical basis, c for Cartesian basis), type of atoms, unknown flag, unknown flag
2. Atom number, number of this type of atoms, number of l types, number of s orbitals, number of p orbitals, ...
3. Atom coordinate
4. H marking the head of a basis, number of alphas, number of contracted basis. Basis are sorted by l
5. Alpha value, contraction coefficients

### MRCI
#### Computation cost
The upper limit of CI expansion on 24 core avx2 processor computer in 2019 is ~ 200,000,000

#### General modification
Tighten MCSCF tolerance in mcscfin before running MRCI:
* energy change < tol(1)=1e-12
* gradient (W) norm < tol(2)=1e-8
* orbital correction (K) norm < tol(3)=1e-8

To run parallel MRCI:
* you may have to manually change the 'ciudg' keyword in control.run to 'pciudg'
* modify nseg in ciudgin according to WORK/ciudg.perf to speed up

Tighten MRCI tolerance in ciudgin before computing gradient:
* RTOLBK=1e-5 or 1e-6
* RTOLCI=1e-5 or 1e-6

Tighten MRCI gradient tolerance in cigrdin:
* orbital resolution < wndtol=1e-9,wnatol=1e-9,wnvtol=1e-9
* final effective densitymatrix < ftol=1e-11
* if using LAPACK solver (no symmetry): solvelpck=1,mdir=0,cdir=0; else Columbus solver (with symmetry): nmiter=200,nvrsmx=200,rtol=1e-10,dtol=1e-10

#### Special skill
* FROOT: NROOT computes the lowest NROOT MRCI roots; FROOT follows the FROOT-th MCSCF root and optimizes it with MRCI (keeping tracking by overlap), then output it along with lower MRCI states. FROOT is better at giving the state you want, but may fail when MRCI has different state ordering from MCSCF

#### Common bug
'maxbl' too small: increase maxbl in cisrtin, usually up to 200000

### Geometry optimization
* NROOT also specifies which surface to optimize on. You can overwrite it with the transition moment input transmomin (i.e. set computing transition moment between m-th and m-th state to optimize on m-th surface)
* By default GDIIS searches for minimum with BFGS. To search for saddle point, replace bfgs at the last line in gdiisin with sadd

### Weird stuff
Orbital
* The number of each kind of orbitals appearing in cidrtin cannot exceed 256, including the frozen orbitals

Interactive input (colinp)
* In CI input, display bug occurs for font size > 12
* In hessian input, perl bug occurs for font size > 8
