## README ##
To start working on a feature:  

1) Clone the repo  
2) Create a local branch and name it after the issue on the tracker (e.g. Issue18)  
3) Work on the issue updating the issue tracker with comments to keep everyone informed of progress  
4) When ready make a commit on your local branch describing the changes in as much detail as possible  
5) Push the changes to a remote version of the repo for all to see  
6) Mark the issue as resolved  
7) Register a pull request so we can discuss the changes before merging into the master branch  
8) Revise or merge as a group into the master  
9) Mark the issue as closed.  

Any questions contact Adrian.

## SETTING UP CODE::BLOCKS WITH MPI LIBRARIES ##
MPI libraries are required in order to compile the code. The steps to installing these (in Ubuntu) and then setting up Codeblocks to include them are given below, as well as the compilation and execution steps.

1) Install MPICH2 using the command: **sudo apt-get install mpich2**  
2) Check it is installed correctly by typing: **mpichversion** (if you get an output then it is installed).  
3) Open your Codeblocks project.  
4) Go to *Settings -> Compiler* and make sure the **GNU GCC Compiler** is selected. Copy the compiler settings and give a name to the new compiler.  
5) Make sure the new compiler is selected and open the **Toolchain Executables** tab. If you kept the default mpich2 installation settings the installation directory should be **/usr**.  
6) Change the C compiler to **mpicc** and the C++ compiler and dynamic linker to **mpicxx**. Press OK.  
7) Go to *Project -> Build Options* and make sure the **Selected Compiler** is set to the one you just created. Note that you have to do this for each individual target you are working on (e.g. Debug, Release). Press OK.  
8) You will likely receive a warning message recommending to rebuild the project as the compiler has been changed. To rebuild go to *Build -> Rebuild*. Note that to compile in serial you can comment out the **BUILD_FOR_MPI** flag in the **definitions.h** header file.  
9) To run the program open up a terminal and navigate to the directory containing the executable. If the program was compiled in serial then the **./LatBo** command will work. However, if it was compiled for parallel execution then use: **mpirun -np N ./LatBo** where N is the number of processes (cores) that you want to use. Note that the number you select for N must match the number of cores you defined during compile time in the **definitions.h** file.

There is probably a way to execute in parallel directly from Codeblocks, when a solution for this is found it will be posted up here. Any questions or issues then contact Joe.

## CODE::BLOCKS RETIRED IN FAVOUR OF ECLIPSE ##
Contact Joe for guidance on setting up.

## SETTING UP VISUAL STUDIO ON WINDOWS ##
Most of the original cost has been written on Windows and compiled using VC++ 2012. Although platform independence is maintained as best as possible there are limitations to maintain compatibility for VC++ 2012 and MSMPI v7 which we used to develop the code originally.
For advice and guidance setting up the code on Visual Studio contact Adrian.

## CONTACTS ##
Adrian Harwood (adrian.harwood@manchester.ac.uk)  
Joe O'Connor (joseph.oconnor@manchester.ac.uk)