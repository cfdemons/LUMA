#pragma once
#include "Marker.h"

// Definition for the BFLMarker class which represents a marker which defines a BFL object //
class BFLMarker :
	public Marker
{

	// Allow BFLbody to access the marker positions etc.
	friend class BFLBody;

public:
	// Default constructor and destructor
	BFLMarker(void);
	~BFLMarker(void);

	// Constructor
	BFLMarker(double x, double y, double z);

protected:


	/*
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/





	/*
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/




};

