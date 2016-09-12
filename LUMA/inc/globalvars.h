/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) 2015, 2016
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * distribution without written consent.
 *
 */

// Header guard
#ifndef GLOBAL_H
#define GLOBAL_H

// Include definitions
#include "definitions.h"

// Global variable references
extern const int c[3][L_nVels];				///< Lattice velocities
extern const double w[L_nVels];				///< Quadrature weights
extern const double cs;						///< Lattice sound speed

#endif
