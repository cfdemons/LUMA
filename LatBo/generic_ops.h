// Header file to reference generic functions
#ifndef HEADER_OPS_GEN
#define HEADER_OPS_GEN

// Forward declarations
std::vector<int> onespace(int min, int max);						// Function: onespace
std::vector<double> linspace(double min, double max, int n);		// Function: linspace
double vecnorm(double vec[2]);										// Function: vecnorm + 4 overloads
double vecnorm(double val1, double val2);
double vecnorm(double val1, double val2, double val3);
double vecnorm(double vec[3]);
double vecnorm(std::vector<double>& vec);
std::vector<int> indmapref(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start); // Function: indmapref
double dotprod(std::vector<double> vec1, std::vector<double> vec2);		// Function: dotprod
std::vector<double> matrix_multiply(std::vector< std::vector<double> >& A, std::vector<double>& x);	// Function: matrix_multiply

#endif