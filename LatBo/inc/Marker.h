#pragma once

#include <vector>

/** Represents a generic marker for a body **/
class Marker
{

public:
	Marker(void)
	{
	};
	~Marker(void)
	{
	};

	// Constructor
	Marker(double x, double y, double z)
	{
		position.push_back(x);
		position.push_back(y);
		position.push_back(z);

	}

	// Protected Members //
protected:
	std::vector<double> position;		// Position vector of marker location in physical units

	// Vector of indices for the LBM nodes considered to be in support of marker
	std::vector<unsigned int> supp_i;
	std::vector<unsigned int> supp_j;
	std::vector<unsigned int> supp_k;


};

