# Columbus MRCI experience

## Computation cost
* The upper limit of CI expansion is ~ 200,000,000

## General modification

* Tighten MCSCF tolerance in mcscfin before running MRCI:
1. energy change < tol(1)=1e-15
2. gradient (W) norm < tol(2)=1e-10
3. orbital correction (K) norm < tol(3)=1e-10
4. start using CSF Hessian when gradient norm < tol(9)=1e-5

* To run parallel CI, you have to manually change the 'ciudg' keyword in control.run to 'pciudg'

* Tighten MRCI tolerance in ciudgin before computing gradient:
1. RTOLCI=1e-5 or 1e-6

* Tighten MRCI gradient tolerance in cigrdin:
1. Orbital resolution < wndtol=1e-9,wnatol=1e-9,wnvtol=1e-9
2. Final effective densitymatrix < ftol=1e-11
3. If using LAPACK solver (no symmetry): solvelpck=1,mdir=0,cdir=0; else Columbus solver (with symmetry): rtol=1e-10,dtol=1e-10,nmiter=200,nvrsmx=200

## Special skill
* FROOT: NROOT computes the lowest NROOT MRCI roots; FROOT follows the FROOT-th MCSCF root and optimizes it with MRCI (keeping tracking by overlap), then output it along with lower MRCI states. FROOT is better at giving the state you want, but may fail when MRCI has different state ordering from MCSCF

## Common bug
* 'maxbl' too small: increase maxbl in cisrtin, usually up to 200000