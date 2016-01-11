#pragma once

// Markers may be of different types depending on the body
template <typename MarkerType>

/** Represents a general body **/
class Body
{

	// Default Constructor / Destructor //
public:
	Body(void)
	{
	};
	~Body(void)
	{
	};

	// Protected Members //

protected:
	double spacing;							// Spacing of the markers in physical units
	std::vector<MarkerType> markers;		// Array of markers which make up the body
	bool closed_surface;					// Flag to specify whether or not it is a closed surface (for output)

	// Methods //
	void addMarker(double x, double y, double z)	// Add marker to the body
	{
		markers.emplace_back(x, y, z);	// Add a new marker object to the array
	}

};

