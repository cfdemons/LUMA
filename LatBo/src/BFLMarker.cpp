#include "..\inc\stdafx.h"
#include "..\inc\BFLMarker.h"
#include "..\inc\BFLBody.h"

// Implementation for the BFLMarker class //


// Default constructor / destructor
BFLMarker::BFLMarker(void)
{
}


BFLMarker::~BFLMarker(void)
{
}


BFLMarker::BFLMarker(double x, double y, double z)
{
	// Add its position
	position.push_back(x);
	position.push_back(y);
	position.push_back(z);

	// Add indices
	std::vector<int> vox = BFLBody::getVoxInd(x,y,z);
	this->supp_i.push_back(vox[0]);
	this->supp_j.push_back(vox[1]);
	this->supp_k.push_back(vox[2]);

}
