/* This file contains all the generic functions for use in the main routine.
In future, we could turn this into its own class with the functions as methods.
*/

#include "stdafx.h"
#include <vector>		// Dynamic array
#include <sstream>		// String stream package
#include "definitions.h"

using namespace std;

// ***************************************************************************************************

// Returns a vector with n uniformly spaced values between min and max
vector<double> linspace(double min, double max, int n)
{
	// Decalre resulting vector
	vector<double> result;

	// Set counter to zero
	int count = 0;
 
	// Loop 
	for (int i = 0; i <= n-2; i++)
	{
		double temp = min + i*(max-min)/(floor( (double)n ) - 1); // Cast n to a double to use floor
		result.insert(result.begin() + count, temp); // Insert element
		count += 1;
	}
 
	// Add last element
	result.insert(result.begin() + count, max);
	
	// Return vector
	return result;
}

// ***************************************************************************************************

// Like linspace but spaces elements by 1
vector<int> onespace(int min, int max)
{
	// Declare resulting array
	vector<int> result;

	// Loop and insert elements
	for (int i = 0; i <= (max-min); i++)
	{
		result.insert( result.begin() + i, min + i ); // Insert element

	}
 
	// Return vector
	return result;
}

// ***************************************************************************************************

// Functions to compute the magnitude of a vector. Overloads allow different input types
// 2D vector with arguments supplied separately
double vecnorm( double val1, double val2 )
{
	double result;
	
	result = sqrt( pow(val1,2) + pow(val2,2) );

	return result;
}

// 3D vector with arguments supplied separately
double vecnorm( double val1, double val2, double val3 )
{
	double result;
	
	result = sqrt( pow(val1,2) + pow(val2,2) + pow(val3,2) );

	return result;
}

// Supplied as a vector
double vecnorm( double vec[] )
{
	double result;
	
#if (dims == 3)

		result = sqrt( pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2) );

#else

		result = sqrt( pow(vec[0],2) + pow(vec[1],2) );

#endif

	return result;
}

// ***************************************************************************************************