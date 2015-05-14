// Header file to reference generic functions
#ifndef HEADER_OPS_GEN
#define HEADER_OPS_GEN

// Forward declarations
std::vector<int> onespace(int min, int max);						// Function: onespace
std::vector<double> linspace(double min, double max, int n);		// Function: linspace
double vecnorm(double vec[2]);										// Function: vecnorm + 3 overloads
double vecnorm(double val1, double val2);
double vecnorm(double val1, double val2, double val3);
double vecnorm(double vec[3]);

#endif