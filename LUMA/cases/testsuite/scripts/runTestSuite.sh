#!/bin/bash

# This script will compile each definition file of the test suite, run it, then check the ouput against the base results
# If you run "this_script.sh > out.dat &" it will run it in the background and print stdout to out.dat file


# Compiler flags and directories (should change this to suit system)
DIR_HDF5=/usr/lib/x86_64-linux-gnu/hdf5/mpich   # Must set this yourself to where the HDF5 library is installed on your system
CC=mpicxx										# Compiler command
CFLAGS="-std=c++0x -O3"							# Compiler flags - don't need -w flag but the warnings will pollute stdout
EXE=LUMA										# Executable
DIR_INC=${DIR_HDF5}/include						# Include directory
DIR_LIB=${DIR_HDF5}/lib							# Library path
LIB=hdf5										# Library


# Set up the variables containing the paths to directories we need
DIR_WORKING=..							# Working directory
DIR_BASE=${DIR_WORKING}/base			# Base directory containing reference results
DIR_DEF=${DIR_WORKING}/defs				# Definition files which will be copied into the LUMA directory
DIR_INPUT=${DIR_WORKING}/input			# Input directory containing files required for running (pointcloud / inlet etc.)
DIR_RES=${DIR_WORKING}/results			# Results directory where the data will be written out
DIR_LUMA=${DIR_WORKING}/../..			# LUMA directory containing the source files to compile


# File for diff comparison
DIFF_FILE=io_lite.Lev0.Reg0.Rnk0.100.dat


# If running with the clean option it will just delete the results folder
while [ ! $# -eq 0 ]
do
	case "$1" in
		--clean | clean | -c)
			rm -rf ${DIR_RES}
			exit
			;;
	esac
	shift
done


# Delete results directory if it already exists and create it again
rm -rf ${DIR_RES}
mkdir ${DIR_RES}


# Get the number of cases to test and print to screen
NCASES=`find ${DIR_DEF} -maxdepth 1 -type f | wc -l`
printf "\n********** LUMA TEST SUITE **********\n"
printf "Beginning test suite -> there are ${NCASES} cases to test\n\n"


# Loop through all cases
#for CASE_DEF_PATH in ${DIR_DEF}/definitions_000.h
for CASE_DEF_PATH in ${DIR_DEF}/*
do

	# Get the case number by stripping the path
	CASE_NUM=${CASE_DEF_PATH%.*}
	CASE_NUM=${CASE_NUM#*_}


	# Inform user we are starting this case (use expansion to treat CASE_NUM as number rather than string)
	printf "Starting case $((CASE_NUM)) of ${NCASES}...\n"


	# Create a directory for this case in the results directory and copy the definitions file there (in case it needs to be inspected) and the input folder
	CASE_RES_PATH=${DIR_RES}/case${CASE_NUM}
	mkdir ${CASE_RES_PATH}
	cp ${CASE_DEF_PATH} ${CASE_RES_PATH}/.
	cp -r ${DIR_INPUT} ${CASE_RES_PATH}/.


	# Copy the definition file into the inc directory in the LUMA source folder
    cp ${CASE_DEF_PATH} ${DIR_LUMA}/inc/definitions.h


	# Get the path to the soon to be created executable
	CASE_EXE=${CASE_RES_PATH}/${EXE}${CASE_NUM}


	# Starting compiling
	printf "Compiling..."


	# Compile the source files and test for success
	if ${CC} ${CFLAGS} -I${DIR_INC} ${DIR_LUMA}/src/* -o ${CASE_EXE} -L${DIR_LIB} -l${LIB} &> ${CASE_RES_PATH}/compile.log; then
		
		# Compilation was a sucess
		printf "success!\n"

		# Run LUMA
		printf "Running LUMA..."
		
		# Temporailiy change working directory to path where executable is (otherwise LUMA writes output to script folder), run LUMA and test for success
		if (cd ${CASE_RES_PATH} && ./${EXE}${CASE_NUM} > /dev/null); then

			# LUMA completed sucessfully
			printf "success!\n"

			# Get the path to the output directory
			DIR_OUT=`find ${CASE_RES_PATH} -maxdepth 1 -type d -name 'output*'`

			# Checking results
			printf "Runnning a diff on the results..."

			# Perform the diff and collect the output
			if diff -q ${DIR_BASE}/case${CASE_NUM}/${DIFF_FILE} ${DIR_OUT}/${DIFF_FILE} > /dev/null; then

				# Check passed
				printf "success!\n"
				printf "Case $((CASE_NUM)) has passed!\n\n"

			else

				# Check failed
				printf "failed (check diff log file)\n"

				# Write out the diff log
				printf "Writing diff log in case folder..."
				diff -y -W 500 ${DIR_BASE}/case${CASE_NUM}/${DIFF_FILE} ${DIR_OUT}/${DIFF_FILE} > ${CASE_RES_PATH}/diff.log

				# Print status
				printf "Case $((CASE_NUM)) failed!\n\n"
			fi
		else

			# LUMA failed during runtime
			printf "failed (check LUMA log file)...skipping case\n\n"
		fi
	else

		# Case did not compile
		echo "failed (check compile log file)...skipping case\n\n"
	fi
done

printf "Finished test suite\n\n"
