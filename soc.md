# How to calculate MRCI spin-orbit coupling
Columbus cannot calculate the spin-orbit coupling by itself; it requires the molcas AMFI integrals, which in turn requires all calculations to be performed based on molcas AO integrals

## Installation
When installing molcas, add `-DLINALG=MKL` to cmake; otherwise molcas would try to build its own linear algebra library

## Molcas basis set
To define customized basis beyond `colinp`, user has to manually modify `molcas.input`, then `infofl` who contains the number of basis in each irreducible representation

`molcas.input` specifies the basis set, which should be defined in `$MOLCAS/basis_library/`, see molcas documentation or basis set exchange for details

## AO ordering
Usually the AOs are ordered atom by atom, but how a single atom's orbitals are ordered? Different programs have different orders

dalton (hermit):
* all s orbitals first, then all p orbitals, then d then f ...
* for p orbitals, 1st-radial `x, y, z` (Cartesian order), then another radial, then another ...
* for d and higher orbitals, 1st-radial `m = -l, ..., 0, ..., l` (spherical order), then another radial, then another ...

molcas (seward):
* all s orbitals first, then all p orbitals, then d then f ...
* for p orbitals, every radial's `px`, every radial's `py`, every radial's `pz` (Cartesian order)
* for d and higher orbitals, every radial's `m = -l, ..., 0, ..., l` (spherical order)

To reuse results based on dalton AO integrals, we need to reorder AO-space vectors (e.g. mocoef) based on the orderings discussed above

Issue: For scf mo, for d orbital `dalton / molcas = sqrt(1 / 12)`, for f orbital `dalton / molcas = sqrt(3 / 20)`; but for mcscf mo there is no such factor?

## MRCI
If multiple states for each multiplicity are desired, then it is usually applaudable to set `update_mode=10` in `ciudgin`. However, in an soc calculation, instead of using user provided `ciudgin`, `$COLUMBUS/runc` would generate `ciudg.drt*` for each multiplicity, which sets `update_mode=1`.

In order to change that, user has to modify `$COLUMBUS/perlscripts/rel_mod.pm`. The original lines are:
```
changekeyword("$ciuin","$ciuin","update_mode",1);
changekeyword("ciudgin","ciudgin","update_mode",1);
```
Which should be replaced with:
```
changekeyword("$ciuin","$ciuin","update_mode",10);
changekeyword("ciudgin","ciudgin","update_mode",10);
```
