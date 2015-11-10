lammps-extensions
=================
Auxiliary code used in lammps and for postprocessing lammps output data formats

* fix_dump_mesh - dumps into OBJ geometry format (angles are used as triangles)
* fix_ave_spatial - modified ave spatial fix which can write into tec data format. If output file has extension *.tec, 
output file format is tec data. It can be opened with TecPlot (probably, need to rename in *.dat) and with Paraview.
* fix_count_atoms - count atoms in a region, uses a custom communicator to be effective
* region_complement - NOT operation on regions
* region_difference - MINUS operation on regions
* molecule_counter - class which simplifies work with molecules in specified group
* gather_containers - collection of routines used to simpify gathering stl containers from different MPI nodes
* atom2plt.sh - script which converts lammps data files (molecular only) into tec format. It can be read by TecPlot.
Example of input data file is cube.atom, output example is cube.plt.
* restart2obj.py - python script which converts a collection of restart files into obj files.
