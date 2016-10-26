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

// stdafx.cpp : source file that includes just the standard includes
// LUMA.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "../inc/stdafx.h"

// Lattice velocities
#if (L_dims == 3) && defined L_USE_KBC_COLLISION

// D3Q27
const int c[3][L_nVels] =
{
	{ 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 0 },
	{ 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 0 },
	{ 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0 }
};

#elseif(L_dims == 3) && !defined L_USE_KBC_COLLISION

// D3Q19
const int c[3][L_nVels] =
{
	{1,	-1,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  0 },
	{0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1,  0,  0,  0,  0,  0 },
	{0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1,  0 }
};

#else

// D2Q9
const int c[3][L_nVels] =
{
	{ 1, -1, 0, 0, 1, -1, 1, -1, 0 },
	{ 0, 0, 1, -1, 1, -1, -1, 1, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

#endif

// Weights for D2Q9, D3Q19 and D3Q27 models
#if (L_dims == 3) && defined L_USE_KBC_COLLISION

const double w[L_nVels] =
{ 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,
1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
8.0 / 27.0 };

#elseif (L_dims == 3) && !defined L_USE_KBC_COLLISION

const double w[nVels] =
{1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
1.0/3.0};

#else

const double w[L_nVels] =
{ 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0 };

#endif

// Lattice sound speed
const double cs = 1.0 / sqrt(3.0);