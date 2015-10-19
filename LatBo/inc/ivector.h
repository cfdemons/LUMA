#pragma once

#include <vector>
/*
	This header defines the ivector class which is like std::vector but has a defined
	operator() to allow automatic flattening of indices before returning a reference of value at
	indexed location.
*/

// Needs to be able to accept different datatypes so template
template <typename GenTyp>
class ivector :	public std::vector<GenTyp> // Define ivector class which inherits from std::vector
{
public:

	// C++ standard does not automatically inherit constructors
	ivector( )
	{
	}

	~ivector( )
	{
	}

	/* Would like to inherit non-default constructor as required but the following
	actually crashes the compiler which is a compiler bug.
	using std::vector :: vector(size_type n, const value_type& val,
                 const allocator_type& alloc = allocator_type());
				 */
	
	// Instead, define the required constructor myself
	ivector(size_t size, GenTyp val) {

		this->resize(size);

		for (size_t i = 0; i < size; i++) {
			this->operator[] (i) = val;
		}

	}
	


	/*	
		:::: USE OF REFERENCES ::::
		Return reference (not the same as a pointer but similar) of the particular element in the 1D array.
		Reference is an alias for an object. Anbything done to the reference is done to the object itself.
		A pointer is an object itself, a reference is only a reference to an actual object.
		Here we delcare a reference to type GenTyp as the return type.
	*/

	// Define overloaded operator() for different array index flattening
	// 4D array index
	inline GenTyp& operator() (size_t i, size_t j, size_t k, size_t v, size_t j_max, size_t k_max, size_t v_max) {

		// Loop over vel then k then j then i
		unsigned int idx = v + (k*v_max) + (j*v_max*k_max) + (i*v_max*k_max*j_max);

		return this->operator[] (idx);

	}

	// 3D array index
	inline GenTyp& operator() (size_t i, size_t j, size_t k, size_t j_max, size_t k_max) {

		// Loop over k then j then i
		unsigned int idx = k + (j*k_max) + (i*k_max*j_max);

		return this->operator[] (idx);

	}

	// 2D array index
	inline GenTyp& operator() (size_t i, size_t j, size_t j_max) {

		// Loop over j then i
		unsigned int idx = j + (i*j_max);

		return this->operator[] (idx);

	}


};