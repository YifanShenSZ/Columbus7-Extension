# Hessian

To calculate Hessian:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Prepare intcfl (Columbus7 internal coordinate definition) and geom (Columbus7 geometry file containing the geometry to calculate Hessian) in your work directory. Run GenLoop.exe (create geom.all), then go to step 3 in SinglePointBatch
3. Go to Collect, then run FiniteDifference.exe (This step creates Hessian, VibrationalFrequency.txt, NormalMode.txt)

After Hessian and normal modes are obtained, you may scan each normal mode:
1. Run NormalModeLoop.exe (create geom.imag and geom.real), then go to step 3 in SinglePointBatch (please rename to geom.all when bash create)
2. Go to Collect. Optionally, run HarmonicityCheck.exe to make sure the scan is valid

Dependency: my Fortran-Library, as written in MyLib of makefile