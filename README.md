## BEFORE YOU BEGIN ##
To start working on a feature:  

1) Clone the repo.  
2) Create a local branch and name it after the issue on the tracker (e.g. Iss18).  
3) Work on the issue updating the issue tracker with comments to keep everyone informed of progress.  
4) When ready make a commit on your local branch update the releases notes in the *./docs/* directory describing the changes in as much detail as possible. Use existing formatting.  
5) Push the changes to a remote version of the repo for all to see.  
6) Mark the issue as resolved.  
7) Register a pull request so we can discuss the changes before merging into the master branch.  
8) Revise or merge as a group into the pre-release branch.  
9) Mark the issue as closed.  

Any questions contact Adrian.

## COMPILING IN LINUX WITH GCC ##
The following steps will guide you through how to compile and run LUMA using the GCC compiler, with the MPICH and HDF5 libraries from a Linux terminal. To set this up in an IDE is a straightforward extension:

1) Install the GCC compiler, MPICH and HDF5 libraries using the command:
```
#!c++
sudo apt-get install gcc mpich libhdf5-mpich-dev
```

2) While it isn't essential to set the ```HDF5_HOME``` environment variable it makes it more convenient. First find the path to where the HDF5 library is installed - for me this is at ```/usr/lib/x86_64-linux-gnu/hdf5/mpich```
3) Set the ```HDF5_HOME``` environment variable:
```
#!c++
export HDF5_HOME=/usr/lib/x86_64-linux-gnu/hdf5/mpich
```
(this can be added to your *.profile* to save you having to do this every time you open a terminal)
4) When compiling use the MPICH wrapper script (with the C++11 option) and make sure the HDF5 libraries are included and linked:
```
#!c++
mpicxx -std=c++0x -O3 -I{HDF5_HOME}/include path_to_LUMA/src/* -o LUMA -L{HDF5_HOME}/lib -lhdf5
```
(the optimisation flag is optional)
5) If everything has gone as planned there will be an executable in your present working directory called **LUMA**. To run LUMA enter the command:
```
#!c++
mpirun -np NPROCS ./LUMA
```
where *NPROCS* is the number of processes you want to run LUMA with (this value must match the number of processes set in the definition file before compiling or LUMA will crash)

Any questions or issues then contact Joe.

## COMPILING ON CSF WITH HDF5 AND MPI ##
Parallel HDF5 libraries are now available on CSF. To build LUMA 1.2+ on CSF you will need to load the compiler libraries, the MPI libraries and the HDF5 libraries (in that order) using

```
#!c++
module load compilers/gcc/4.9.0
module load mpi/gcc/openmpi/1.6-ib
module load apps/gcc/hdf5/1.8.16-mpi
```

and then LUMA can be built using the command

```
#!c++
mpiCC -w -std=c++11 -I$HDF5DIR/include  src/*.cpp  -o LUMA  -L$HDF5DIR/lib -lhdf5 –llapack
```

Any questions contact Adrian.

## SETTING UP CODE::BLOCKS WITH MPI LIBRARIES ##
MPI libraries are required in order to compile the code. The steps to installing these (in Ubuntu) and then setting up Codeblocks to include them are given below, as well as the compilation and execution steps.

1) Install MPICH2 using the command:
```
#!c++
sudo apt-get install mpich2
```  
2) Check it is installed correctly by typing:
```
#!c++
mpichversion
```
(if you get an output then it is installed).  
3) Open your Codeblocks project.  
4) Go to *Settings -> Compiler* and make sure the **GNU GCC Compiler** is selected. Copy the compiler settings and give a name to the new compiler.  
5) Make sure the new compiler is selected and open the **Toolchain Executables** tab. If you kept the default mpich2 installation settings the installation directory should be **/usr**.  
6) Change the C compiler to **mpicc** and the C++ compiler and dynamic linker to **mpicxx**. Press OK.  
7) Go to *Project -> Build Options* and make sure the **Selected Compiler** is set to the one you just created. Note that you have to do this for each individual target you are working on (e.g. Debug, Release). Press OK.  
8) You will likely receive a warning message recommending to rebuild the project as the compiler has been changed. To rebuild go to *Build -> Rebuild*. Note that to compile in serial you can comment out the ```BUILD_FOR_MPI``` flag in the **definitions.h** header file.  
9) To run the program open up a terminal and navigate to the directory containing the executable. If the program was compiled in serial then the ```./luma``` command will work. However, if it was compiled for parallel execution then use:
```
#!c++
mpirun -np N ./luma
```
where N is the number of processes (cores) that you want to use. Note that the number you select for N must match the number of cores you defined during compile time in the **definitions.h** file.

There is probably a way to execute in parallel directly from Codeblocks, when a solution for this is found it will be posted up here. Any questions or issues then contact Joe.

## CODE::BLOCKS RETIRED IN FAVOUR OF ECLIPSE ##
Contact Joe for guidance on setting up.

## SETTING UP VISUAL STUDIO ON WINDOWS ##
Originally, the code was written on Windows and compiled with VC++ 2012. Although platform independence is maintained as best as possible there are limitations to maintain compatibility for the somewhat "wild" VC++ 2012 and MSMPI v7. The software has also been successfully compiled using VC++ 2013 which is the current development platform on Windows. In addition, users will also need Parallel HDF 1.8.17 installed and for the merge tool compilation will need VTK 7.0.
For advice and guidance setting up the code and compiling on Windows with Visual Studio contact Adrian.

## CONTACTS ##
Adrian Harwood (adrian.harwood@manchester.ac.uk)  
Joe O'Connor (joseph.oconnor@manchester.ac.uk)