/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

// stdafx.cpp : source file that includes just the standard includes
// LUMA.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "../inc/stdafx.h"

// Lattice velocities
#if (L_DIMS == 3) && defined L_USE_KBC_COLLISION

// D3Q27
const int c[3][L_NUM_VELS] =
{
	{ 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 0 },
	{ 0, 0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 0 },
	{ 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0 }
};

const int c_opt[L_NUM_VELS][3] =
{
	{ 1, 0, 0 },
	{ -1, 0, 0 },
	{ 0, 1,	0 },
	{ 0, -1, 0 },
	{ 0, 0,	1 },
	{ 0, 0,	-1 },
	{ 0, 1,	1 },
	{ 0, -1, -1 },
	{ 0, 1, -1 },
	{ 0, -1, 1 },
	{ 1, 0,	1 },
	{ -1, 0, -1 },
	{ 1, 0,	-1 },
	{ -1, 0, 1 },
	{ 1, 1,	0 },
	{ -1, -1, 0 },
	{ 1, -1, 0 },
	{ -1, 1, 0 },
	{ 1, 1,	1 },
	{ -1, -1, -1 },
	{ -1, -1, 1, },
	{ 1, 1,	-1 },
	{ -1, 1, 1 },
	{ 1, -1, -1 },
	{ 1, -1, 1 },
	{ -1, 1, -1 },
	{ 0, 0,	0 }
};

#elif(L_DIMS == 3) && !defined L_USE_KBC_COLLISION

// D3Q19
const int c[3][L_NUM_VELS] =
{
	{1,	-1,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  0 },
	{0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1,  0,  0,  0,  0,  0 },
	{0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1,  0 }
};

const int c_opt[L_NUM_VELS][3] =
{
	{ 1, 0, 0 },
	{ -1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, -1, 0 },
	{ 0, 0, 1 },
	{ 0, 0, -1 },
	{ 1, 1, 0 },
	{ -1, -1, 0 },
	{ 1, -1, 0 },
	{ -1, 1, 0 },
	{ 0, 1, 1 },
	{ 0, -1, -1 },
	{ 0, 1, -1 },
	{ 0, -1, 1 },
	{ 1, 0, 1 },
	{ -1, 0, -1 },
	{ -1, 0, 1 },
	{ 1, 0, -1 },
	{ 0, 0, 0, }
};

#else

// D2Q9
const int c[3][L_NUM_VELS] =
{
	{ 1, -1, 0, 0, 1, -1, 1, -1, 0 },
	{ 0, 0, 1, -1, 1, -1, -1, 1, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

const int c_opt[L_NUM_VELS][3] =
{
	{ 1, 0, 0 },
	{ -1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, -1, 0 },
	{ 1, 1, 0 },
	{ -1, -1, 0 },
	{ 1, -1, 0 },
	{ -1, 1, 0 },
	{ 0, 0, 0, }
};

#endif

// Weights for D2Q9, D3Q19 and D3Q27 models
#if (L_DIMS == 3) && defined L_USE_KBC_COLLISION

const double w[L_NUM_VELS] =
{ 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0, 2.0 / 27.0,
1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0, 1.0 / 54.0,
1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0, 1.0 / 216.0,
8.0 / 27.0 };

#elif (L_DIMS == 3) && !defined L_USE_KBC_COLLISION

const double w[L_NUM_VELS] =
{1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
1.0/3.0};

#else

const double w[L_NUM_VELS] =
{ 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0 };

#endif

// Lattice sound speed
const double cs = 1.0 / sqrt(3.0);