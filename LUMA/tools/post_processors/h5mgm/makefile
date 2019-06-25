# Makefile for merge tool (will probably need to edit compiler commands and paths for specific machine and install)

# Compiler command
CC=mpicxx -O3 -std=c++0x -w

# VTK Version
VTK_VER=6.2

# HDF5 paths
DIR_HDF5=${HDF5_HOME}
INC_HDF5=-I${DIR_HDF5}/include
LIBPATH_HDF5=-L${DIR_HDF5}/lib
LIB_HDF5=-lhdf5

# VTK paths
INC_VTK=-I/usr/include/vtk-$(VTK_VER)
LIBPATH_VTK=-L/usr/lib/x86_64-linux-gnu
LIB_VTK=-lvtkFiltersCore-$(VTK_VER) -lvtkIOLegacy-$(VTK_VER) -lvtkIOCore-$(VTK_VER) -lvtkCommonExecutionModel-$(VTK_VER) -lvtkCommonDataModel-$(VTK_VER) -lvtkCommonCore-$(VTK_VER) -lvtkIOXML-$(VTK_VER)


# Location of source, header and object files
SDIR=./src
HDIR=.
ODIR=.
TDIR=$(SDIR)/ThirdParty


# List of header files
DEPS = h5mgm.h
TDEPS = $(TDIR)/vtkCleanUnstructuredGrid.h


# List of object files
OBJ = h5mgm.o
TOBJ = $(TDIR)/vtkCleanUnstructuredGrid.o


# Compile the source files into object files
$(TDIR)/%.o: $(TDIR)/%.cpp $(TEPS)
	$(CC) -c -o $@ $< $(INC_HDF5) $(INC_VTK)


# Compile the source files into object files
$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(INC_HDF5) $(INC_VTK)


# Link object files to get executable
h5mgm: $(OBJ) $(TOBJ)
	$(CC) -o $@ $^ $(LIBPATH_HDF5) $(LIB_HDF5) $(LIBPATH_VTK) $(LIB_VTK)


# Clean the project
.PHONY: clean

# Clean up the directory
clean:
	rm -rf *.o $(TDIR)/*.o h5mgm
