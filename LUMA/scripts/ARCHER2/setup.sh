#!/bin/bash

module load cray-hdf5-parallel

export MPICXX="cc"
export LIB="-lhdf5 -lstdc++"
