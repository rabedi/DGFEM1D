#ifndef ONE_DIMENSIONAL_ELEMENT__H
#define ONE_DIMENSIONAL_ELEMENT__H

#include "LAFuncs.h"
#include "globalTypesClasses.h"

// this is the parent element from xi = -1 to 1 and no real geometry or material properties are considered here
// N is shape function
// Bxi = dN / dxi
class OneDimensionalParentElement
{
public:
	// element order
	unsigned int polyOrder;
	// whether we want the mass-lumpting option or not
	bool lumpMass;

	void Initialize(unsigned int polyOrderIn, bool lumpMassIn);

	/// computed
	unsigned int ndof;
	MATRIX mpe, kpe; // me = int_xi=-1 -> 1 N^t N d xi,   ke = int_xi=-1 -> 1 Bxi^t Bxi d xi
};


class OneDimensionalElement
{
public:
	ElementProperties elementProps;
	MATRIX me, ce, ke;
};

#endif