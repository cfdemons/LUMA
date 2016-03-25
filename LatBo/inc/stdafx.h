// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// Frequently used headers (speeds up compilation in VS if put in the pre-compiled header module)
#include <algorithm>
#include <cmath>
#include <vector>

// Check OS is Windows or not
#ifdef _WIN32

#include "targetver.h"
#define NOMINMAX	// Stop Windows.h redefining min/max
#include <Windows.h>
#include <tchar.h>

#else // Compiling with gcc through Code::Blocks on Linux

#include <stdlib.h> // Includes exit() function
#include <cstring>

#endif

#include <stdio.h>

// Grid utilities class definition (available to all parts of code)
#include "../inc/GridUtils.h"

// Function: is_nan
template <typename NumType>
inline static bool is_nan(NumType n) {

	// Try test so it is platform independent
	if (
#ifdef _WIN32
		_isnan(n)
#else
		isnan(n)
#endif
	) return true;
				
	else return false;

};