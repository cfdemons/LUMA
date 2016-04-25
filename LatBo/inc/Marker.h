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

	// Members //
public:
	std::vector<double> position;		// Position vector of marker location in physical units

	/* Vector of indices for lattice sites considered to be in support of the marker:
	 * In IBM these are all the lattice support sites for spreading and interpolating;
	 * In BFL these are the voxel indices in which the BFL marker resides;
	 */
	std::vector<int> supp_i;
	std::vector<int> supp_j;
	std::vector<int> supp_k;

	// Array of indices indicating on which rank the given support point resides
	std::vector<int> support_rank;


};

