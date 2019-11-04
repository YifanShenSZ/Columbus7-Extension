# Hessian

To calculate Hessian:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Entre a work directory
3. Prepare intcfl (Columbus7 internal coordinate definition file) and geom (Columbus7 geometry file) in current directory. Run GenLoop.exe (create geom.all), then go to step 5 in SinglePointBatch
4. Go to Collect, then run FiniteDifference.exe

After Hessian and normal modes are obtained, you may scan each normal mode:
1. Run NormalModeLoop.exe (create geom.imag and geom.real), then go to step 5 in SinglePointBatch
2. Go to Collect

Dependency: my Fortran-Library, as written in MyLib of makefile
