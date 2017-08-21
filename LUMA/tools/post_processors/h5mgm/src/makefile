# Makefile for merge tool (will probably need to edit compiler commands and paths for specific machine and install)

# Compiler command
CC=mpicxx -O3 -std=c++0x -w

# HDF5 paths
DIR_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/mpich
INC_HDF5=-I${DIR_HDF5}/include
LIBPATH_HDF5=-L${DIR_HDF5}/lib
LIB_HDF5=-lhdf5

# VTK paths
INC_VTK=-I/usr/include/vtk-6.2
LIBPATH_VTK=-L/usr/lib/x86_64-linux-gnu
LIB_VTK=-lvtkFiltersCore-6.2 -lvtkIOLegacy-6.2 -lvtkIOCore-6.2 -lvtkCommonExecutionModel-6.2 -lvtkCommonDataModel-6.2 -lvtkCommonCore-6.2


# Location of source, header and object files
SDIR=.
HDIR=.
ODIR=.


# List of header files
DEPS = h5mgm.h


# List of object files
OBJ = h5mgm.o


# Compile the source files into object files
$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(INC_HDF5) $(INC_VTK)


# Link object files to get executable
h5mgm: $(OBJ)
	$(CC) -o $@ $^ $(LIBPATH_HDF5) $(LIB_HDF5) $(LIBPATH_VTK) $(LIB_VTK)


# Clean the project
.PHONY: clean

# Clean up the directory
clean:
	rm -rf *.o h5mgm
