
INCLUDES = -I/usr/local/include -I/Users/kirill/bin/hdf5/include -I/Users/kirill/Documents/workspace/lammps-rbc/lammps-rbc/lib/levelset-light/include
LIBS_PATH = -L/Users/kirill/bin/hdf5/lib
LIBS = -lhdf5 -lhdf5_cpp

txt2hdf5:	main.cpp
	g++-6  main.cpp -g3 -O0 -fopenmp -std=c++0x -DUSE_HDF5  $(INCLUDES) $(LIBS_PATH) $(LIBS) -o txt2hdf5
	
clean:
	rm txt2hdf5
	
.PHONY = clean
