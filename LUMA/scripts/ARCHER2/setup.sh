#!/bin/bash

module load cray-hdf5-parallel

export CC="cc"
export LIB="-lhdf5 -lstdc++"
