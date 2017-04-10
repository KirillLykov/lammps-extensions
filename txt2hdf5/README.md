txt2hdf5
=================

This is a command line tool to convert output of LAMMPS ave/spatial output text file to hdf5.
This allows to analyse the data in Paraview (open accompanying xmf instead of h5 itself).
In order to compile this code, one should have boost as well as levelset-light (https://github.com/KirillLykov/levelset-light) 
which is header-only library to work with different formats.

Example:
In case if you simulate a flow problem in LAMMPS, write in your script:

compute cc2 all chunk/atom bin/3d x center 0.5 y center 0.5 z center 0.5
fix profile2 all ave/chunk 10 10 100 cc2 vx vy vz norm all file vel3d.txt

As the result of execution, you will get vel3d.txt. This file can contain several frames (niterations/100 in our case), so 
if we want to examine the firts frame we specify:

./txt2hdf5 vel3d.txt outFileName 1 0.5 0.5 0.5

Note, that we also specified the bin size (0.5 for all three dimentions).

