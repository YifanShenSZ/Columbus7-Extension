# The bugs found in Columbus7
Columbus7 originates in 1980s and preserves many historical codes. Besides, the user community of Columbus7 is always small. As a result, there are compatibility issues and undiscovered bugs

## Installation
ifort 2018.3 with avx2 leads to bug in `mcscf.x`

## Interactive input (colinp)
In CI input, display bug occurs for font size > 12

In hessian input, perl bug occurs for font size > 8

## MCSCF
Freezing mcscf orbitals causes erroneous Hamiltonian matrix elements

Energies may differ by 1e-6 Hatree from calculations with or without symmetry (although the molecular orbitals look the same)

The molecular orbitals converged by higher symmetry may not converge at the 1st iteration in lower symmetry

## MRCI
Switching orbital ordering in `cidrtin` would give different number of CSFs, although the energies seem to be the same

Increase `maxbl` in `cisrtin` (usually to 200000), otherwise "maxbl too small" exception may be raised

Frozen virtual is not supported in `cigrd.x`

## Gradient
If specifying `cigrad` in `control.run`, then `runc` assumes geometry optimization job, so `gdiisin` is required

If specifying `cigrad` in `control.run`, then `cigrdin` cannot contain `nadcalc`

If specifying `nadcoupl` in `control.run`, then sometimes calculating only energy gradients without nonadiabatic couplings would make `cigrd.x` raise inconsistent-energy error (probably buggy density matrix)

## Internal coordinate
Out of plane a-b-c-d actually has d as centre, not as claimed in https://www.univie.ac.at/columbus/docs_COL70/fulldoc/intc.txt
