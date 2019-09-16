# Columbus tricks

* Installation
1. You have to copy GA static libraries to Columbus source code directory during installation
2. Better use openmpi rather than intel mpi

* Job control
1. Single point calculation: Never use NROOT to specified the number of states to compute, use FROOT instead (i.e. set both root and follow in colinp)
2. Geometry optimization: NROOT specifies which surface to optimize on. You can overwrite it with the transition moment input (i.e. set computing transition moment between n-th and n-th state)
2. Minimum energy crossing point search: use modified polyhess on Alex2

* Computation cost
1. CI expansion should not exceed 100,000,000

* Historical issue
1. The number of each kind of orbitals appearing in cidrtin cannot exceed 256, including the frozen orbitals! (It's OK to code with unsigned short, but why you enforce it in every aspect??)

* Error message
1. When you see "drt size error" in ciudg execution, it is arised from slice > 6 or strtdrt + forbyt(bufsiz) - 1 > forbyt(drt%drtsiz)