// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// Check OS is Windows or not
#ifdef _WIN32

#include "targetver.h"
#include <tchar.h>

#else // Compiling with gcc through Code::Blocks on Linux

#include <stdlib.h> // Includes exit() function

#endif

#include <stdio.h>
