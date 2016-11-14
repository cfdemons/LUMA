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

# Master log file
LOG_FILE=testsuite.log

# File for diff comparison
DIFF_FILE=io_lite.Lev0.Reg0.Rnk0.100.dat

# If running with the clean option it will just delete the results folder
OPTION="none"
while [ ! $# -eq 0 ]
do
	case "$1" in
		--clean)
			rm -rf ${DIR_RES}
			rm -f ${DIR_WORKING}/${LOG_FILE}
			exit
			;;
		--rebase)
			OPTION="rebase"
			rm -rf ${DIR_BASE}
			mkdir ${DIR_BASE}
			;;
	esac
	shift
done

# Delete results directory if it already exists and create it again
rm -rf ${DIR_RES}
mkdir ${DIR_RES}
rm -f ${DIR_WORKING}/${LOG_FILE}

# Get the number of cases to test and print to screen
NCASES=`find ${DIR_DEF} -maxdepth 1 -type f | wc -l`
printf "\n********** LUMA TEST SUITE **********\n"

# Print out starting message
if [ ${OPTION}  == "rebase" ]; then
	printf "Beginning rebasing -> there are ${NCASES} cases to rebase\n\n"
elif [ ${OPTION} == "none" ]; then
	printf "Beginning test suite -> there are ${NCASES} cases to test\n\n"
fi

# Print header for master log file
DATE=`date +%Y-%m-%d:%H:%M:%S`
printf "\n********** LUMA TEST SUITE - ${DATE} **********\n\n" > ${DIR_WORKING}/${LOG_FILE}

# Loop through all cases
for CASE_DEF_PATH in ${DIR_DEF}/*
do

	# Get the case number by stripping the path
	CASE_NUM=${CASE_DEF_PATH%.*}
	CASE_NUM=${CASE_NUM#*_}
	CASE_NUM_INT=$((10#${CASE_NUM}))

	# Inform user we are starting this case (use expansion to treat CASE_NUM as number rather than string)
	printf "Starting case ${CASE_NUM_INT} of $((NCASES-1))...\n"

	# Create a directory for this case in the results directory and copy the definitions file there (in case it needs to be inspected)
	CASE_RES_PATH=${DIR_RES}/case${CASE_NUM}
	mkdir ${CASE_RES_PATH}
	cp ${CASE_DEF_PATH} ${CASE_RES_PATH}/.

	# Copy the proper inputs to the input directory
	# IB point cloud cases
	if [ ${CASE_NUM_INT} -eq 0 ] || [ ${CASE_NUM_INT} -eq 1 ]; then
		mkdir -p ${CASE_RES_PATH}/input
		cp ${DIR_INPUT}/ibb_input2D.in ${CASE_RES_PATH}/input/ibb_input.in
	fi
	if [ ${CASE_NUM_INT} -eq 30 ] || [ ${CASE_NUM_INT} -eq 31 ]; then
		mkdir -p ${CASE_RES_PATH}/input
		cp ${DIR_INPUT}/ibb_input3D.in ${CASE_RES_PATH}/input/ibb_input.in
	fi

	# BB point cloud cases
    if [ ${CASE_NUM_INT} -eq 3 ] || [ ${CASE_NUM_INT} -eq 6 ]; then
		mkdir -p ${CASE_RES_PATH}/input
		cp ${DIR_INPUT}/bbb_input2D.in ${CASE_RES_PATH}/input/bbb_input.in
	fi
	# BB point cloud cases
	if [ ${CASE_NUM_INT} -eq 33 ] || [ ${CASE_NUM_INT} -eq 36 ]; then
		mkdir -p ${CASE_RES_PATH}/input
		cp ${DIR_INPUT}/bbb_input3D.in ${CASE_RES_PATH}/input/bbb_input.in
	fi

	# BFL point cloud cases
    if [ ${CASE_NUM_INT} -eq 4 ] || [ ${CASE_NUM_INT} -eq 7 ]; then
		mkdir -p ${CASE_RES_PATH}/input
		cp ${DIR_INPUT}/bfl_input2D.in ${CASE_RES_PATH}/input/bfl_input.in
	fi
	if [ ${CASE_NUM_INT} -eq 34 ] || [ ${CASE_NUM_INT} -eq 37 ]; then
		mkdir -p ${CASE_RES_PATH}/input
		cp ${DIR_INPUT}/bfl_input3D.in ${CASE_RES_PATH}/input/bfl_input.in
	fi

	# Read in inlet profile
	if [ ${CASE_NUM_INT} -eq 3 ] || [ ${CASE_NUM_INT} -eq 6 ] || [ ${CASE_NUM_INT} -eq 33 ]; then
		mkdir -p ${CASE_RES_PATH}/input
		cp ${DIR_INPUT}/inlet_profile2D.in ${CASE_RES_PATH}/input/inlet_profile.in
	fi

	# Copy restart data from previous case
    if [ ${CASE_NUM_INT} -eq 1 ] || [ ${CASE_NUM_INT} -eq 9 ] || [ ${CASE_NUM_INT} -eq 31 ]; then
	    mkdir -p ${CASE_RES_PATH}/input
		LAST_CASE=$((${CASE_NUM_INT}-1))			
		LAST_CASE=$(printf %03d ${LAST_CASE})
		if [ ${OPTION}  == "rebase" ]; then
			LAST_DIR_OUT=${DIR_BASE}/case${LAST_CASE}
		elif [ ${OPTION} == "none" ]; then
			LAST_DIR_OUT=`find ${DIR_RES}/case${LAST_CASE} -maxdepth 1 -type d -name 'output*'`
		fi
		cp ${LAST_DIR_OUT}/restart* ${CASE_RES_PATH}/input/.
	fi

	# Number of processes to run on
	if [ ${CASE_NUM_INT} -eq 0 ] || [ ${CASE_NUM_INT} -eq 1 ] || [ ${CASE_NUM_INT} -eq 2 ] || [ ${CASE_NUM_INT} -eq 3 ] || [ ${CASE_NUM_INT} -eq 4 ] || [ ${CASE_NUM_INT} -eq 5 ] || [ ${CASE_NUM_INT} -eq 30 ] || [ ${CASE_NUM_INT} -eq 31 ] || [ ${CASE_NUM_INT} -eq 32 ] || [ ${CASE_NUM_INT} -eq 33 ] || [ ${CASE_NUM_INT} -eq 34 ] || [ ${CASE_NUM_INT} -eq 35 ]; then
		NPROCS=1
	elif [ ${CASE_NUM_INT} -eq 6 ] || [ ${CASE_NUM_INT} -eq 7 ] || [ ${CASE_NUM_INT} -eq 8 ] || [ ${CASE_NUM_INT} -eq 9 ]; then
		NPROCS=8
	fi

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
		if (cd ${CASE_RES_PATH} && mpirun -np ${NPROCS} ./${EXE}${CASE_NUM} > LUMA.log); then

			# LUMA completed sucessfully
			printf "success!\n"

			# Get the path to the output directory
			DIR_OUT=`find ${CASE_RES_PATH} -maxdepth 1 -type d -name 'output*'`

			# If rebasing then move this to the base directory
			if [ ${OPTION} == "rebase" ]; then
				mv ${DIR_OUT} ${DIR_BASE}/case${CASE_NUM}
				mv ${CASE_RES_PATH}/* ${DIR_BASE}/case${CASE_NUM}
				printf "Case ${CASE_NUM_INT} has been rebased successfully!\n\n"
				printf "CASE ${CASE_NUM_INT} -> REBASED SUCCESSFULLY\n" >> ${DIR_WORKING}/${LOG_FILE}
			
			# If not rebasing then perform the diff on the results			
			elif [ ${OPTION} == "none" ]; then

				# Checking results
				printf "Runnning a diff on the results..."

				# Perform the diff and collect the output
				if diff -q ${DIR_BASE}/case${CASE_NUM}/${DIFF_FILE} ${DIR_OUT}/${DIFF_FILE} > /dev/null; then

					# Check passed
					printf "success!\n"
					printf "Case ${CASE_NUM_INT} has passed!\n\n"

					# Write to master log file
					printf "CASE ${CASE_NUM_INT} -> PASS\n" >> ${DIR_WORKING}/${LOG_FILE}

				else

					# Check failed
					printf "failed (check diff log file)\n"

					# Write out the diff log
					printf "Writing diff log in case folder..."
					diff -y -W 500 ${DIR_BASE}/case${CASE_NUM}/${DIFF_FILE} ${DIR_OUT}/${DIFF_FILE} > ${CASE_RES_PATH}/diff.log

					# Print status
					printf "Case ${CASE_NUM_INT} failed!\n\n"

					# Write to master log file
					printf "CASE ${CASE_NUM_INT} -> FAILED ON DIFF\n" >> ${DIR_WORKING}/${LOG_FILE}
				fi
			fi
		else

			# LUMA failed during runtime
			printf "failed (check LUMA log file)...skipping case\n\n"

			# Write to master log file
			printf "CASE ${CASE_NUM_INT} -> FAILED ON RUN\n" >> ${DIR_WORKING}/${LOG_FILE}

		fi
	else

		# Case did not compile
		echo "failed (check compile log file)...skipping case\n\n"

		# Write to master log file
		printf "CASE ${CASE_NUM_INT} -> FAILED ON COMPILE\n" >> ${DIR_WORKING}/${LOG_FILE}

	fi
done

# Write out finished message
if [ ${OPTION}  == "rebase" ]; then
	rm -rf ${DIR_RES}
	printf "Finished rebasing test suite\n\n"
elif [ ${OPTION} == "none" ]; then
	printf "Finished test suite\n\n"
fi
