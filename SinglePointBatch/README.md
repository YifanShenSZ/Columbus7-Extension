# Single point batch

To run a batch of Columbus7 single point calculation:
1. Prepare a directory with input template, i.e. everything required for a single point calculation except geom
2. Modify GenColumbusInput.f90 to generate a batch of desired geometry, make, make clean, then ./GenColumbusInput.exe (This step creates geom.all)
3. Modify the user input section in create (for path of input template directory, and NAtoms), bash create (This step creates directories naming from 0 to NGeoms - 1, which are single point calculation job directories)
4. Modify RunColumbus to fit your queue, then bash submit $1, where $1 = NGeoms (This step submits the batch jobs. If some jobs fail, bash submit $1 again)