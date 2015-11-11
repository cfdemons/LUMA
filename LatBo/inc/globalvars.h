// Header guard
#ifndef GLOBAL_H
#define GLOBAL_H

// Include definitions
#include "definitions.h"

// Global variable references
extern const int c[3][nVels];				// Lattice velocities
extern const double w[nVels];				// Weights
extern const double cs;						// Lattice sound speed for lattice
extern time_t curr_time;                    // Current time returned from time_t = time(NULL)
extern std::string timeout_str;             // The output directory path with the formatted date/time stamp

#ifdef USE_MRT
extern const int mMRT[nVels][nVels];		// MRT transformation matrix
extern const double mInvMRT[nVels][nVels];	// Inverse of MRT transformation matrix
#endif

#endif
