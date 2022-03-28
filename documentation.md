# Complement to Columbus7 documentation
Here are the solutions to the issues I have encounterd, which are out of the coverage of the official documentation

## Installation
Steps
1. Install a *Global Arrays*
2. Enter a desired directory, unzip the installation package
3. cp (path to *Global Arrays*)/lib/*.a Columbus
4. ./install.automatic cpan standard grad parallel
5. export COLUMBUS=(current directory)/Columbus

Solutions to some issues
* Global Arrays (subsequently Columbus7) does not support mpi 3.0 standard, e.g. openmpi 4 fails but 3 works
* Installation may fail with intel mpi. For now openmpi always works
* Replace certain files with their counterparts in Modification

## Basis set
Linear dependency will not be excluded automatically, so user has to take care of that

To define customized basis beyond `colinp`, user has to manually modify `daltaoin`, then `infofl` who contains the number of basis in each irreducible representation

`daltaoin` specifies the GTOs and the contractions, whose format is:
1. (s for spherical basis, c for Cartesian basis), type of atoms, unknown flag, unknown flag
2. Atom number, number of this type of atoms, number of l types, number of s orbitals, number of p orbitals, ...
3. Atom coordinate
4. H marking the head of a basis, number of alphas, number of contracted basis. Basis are sorted by l
5. Alpha value, contraction coefficients

## MCSCF
MCSCF is a nonlinear optimization problem, so may easily be trapped into local minimum, e.g.:
* Molecular orbital may lose symmetry if not explicitly using symmetry (although this may also arise from intruder states)
* Even the molecular orbitals have correct symmetry, the calculation without symmetry may give a higher energy

Tighten MCSCF tolerance in mcscfin before running MRCI, recommend:
* gradient (W) norm < tol(2)=1e-8
* orbital correction (K) norm < tol(3)=1e-8

## MRCI
The upper limit of CI expansion on a 24 CPU 128 GB memory computer in 2020 is ~ 200,000,000

To run parallel MRCI:
* you may have to manually change the 'ciudg' keyword in control.run to 'pciudg'
* modify nseg in ciudgin according to WORK/ciudg.perf to speed up

Tighten MRCI tolerance in ciudgin before computing gradient, recommend:
* RTOLCI = 1e-5

Tighten MRCI gradient tolerance in cigrdin:
* orbital resolution < wndtol = 1e-9, wnatol = 1e-9, wnvtol = 1e-9
* final effective densitymatrix < ftol = 1e-11
* if using LAPACK solver (no symmetry): solvelpck = 1, mdir = 0, cdir = 0; else Columbus solver (with symmetry): nmiter = 200, nvrsmx = 200, rtol = 1e-10, dtol = 1e-10

## Geometry optimization
`gdiis.x` can be used to search for minimum, saddle point, minimum energy crossing by setting keyword `bfgs`, `sadd`, `coni`

The initial Hessian is inputed from:
1. `WORK/hessianinv`, for `bfgs`
2. `hessian`, for `bfgs` and `sadd`
3. `intcfl`, for `bfgs` and `sadd` and `coni`

Only `bfgs` updates Hessian by BFGS method, others keep using initial Hessian

To continue `gdiis.x` iteration, copy `LISTINGS/gdiisfl` to `WORK` since it stores the historical geometires to form the iterative subspace

To continue `coni`, additionally copy `WORK/zetafl` to `WORK` since it stores the Lagrangian multipliers

### Saddle point
RGF is claimed to be a global search for saddle point, but no one has ever tried that since `gdiis.x` works well

### Minimum energy crossing
`polyhess.x` is Yarkony's method. It requires 2 steps:
1. To find a degenerate seam
2. To minimize energy on the seam

To continue `polyhess.x` iteration, copy `WORK/h-pieces` and `WORK/continuity` to `WORK` since they store the BFGS Hessian and the Lagrangian multipliers

`polyhesin` to approach degenerate seam:
* &NACINT{maxit=200,newton=1,iheseq1=1,ihess=0,ipflg=3,accel=1,scale=0.1,kscale=0}/end
* This is actually steepest descent (always use unit Hessian), with learning rate = `scale`
* So there is no need to continue with WORK/h-pieces
* Cris recommends not to continue with WORK/continuity as well

`polyhesin` to start minimizing energy on the seam (building Hessian):
* &NACINT{maxit=200,newton=1,iheseq1=-1,methodn=99*-1,ihess=0,ipflg=3,accel=1,scale=0.1,kscale=0,}/end
* The difference to Newton iteration is learning rate = `scale` < 1

`polyhesin` to end search:
* &NACINT{maxit=200,newton=1,iheseq1=-1,methodn=99*-1,ihess=0,ipflg=3,accel=1,scale=1.0,kscale=2,}/end
* This is Newton iteration

## Flags that do not have to be tuned
Numerical accuracy:
* integral accuracy in `daltaoin`
* AO-MO transformation accuracy in `tranin`

Numerical stability:
* `FRCSUB` in `ciudgin` (if you want to tune it, check AO linear dependency first)
