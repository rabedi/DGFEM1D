#include "OneDimensionalElement.h"

void OneDimensionalParentElement::Initialize(unsigned int polyOrderIn, bool lumpMassIn)
{
	lumpMass = lumpMassIn;
	polyOrder = polyOrderIn;
	ndof = polyOrder + 1;
	kpe.resize(ndof, ndof);
	mpe.resize(ndof, ndof);

	if (polyOrder == 0)
	{
		kpe = 0.0;
		mpe[0][0] = 2.0;
	}
	else if (polyOrder == 1)
	{
		if (!lumpMass)
		{
			double fact = 1.0 / 3.0;
			mpe[0][0] = mpe[1][1] = 2.0 * fact;
			mpe[0][1] = fact;
			mpe[1][0] = fact;
		}
		else
		{
			mpe[0][0] = mpe[1][1] = 1.0;
			mpe[0][1] = 0.0;
			mpe[1][0] = 0.0;
		}
		// stiffness
		kpe[0][0] = kpe[1][1] = 0.5;
		kpe[0][1] = -0.5;
		kpe[1][0] = -0.5;
	}
	else if (polyOrder == 2)
	{
		// check these mass matrices and stiffness
		if (!lumpMass)
		{
			double fact = 1.0 / 60.0;
			mpe[0][0] = mpe[1][1] = 16.0 * fact;
			mpe[0][1] = mpe[1][0] = 4.0 * fact;
			mpe[2][2] = 32.0 * fact;
			mpe[0][2] = mpe[2][0] = mpe[1][2] = mpe[2][1] = 8.0 * fact;
		}
		else
		{
			double fact = 1.0 / 3.0;
			mpe[0][0] = mpe[1][1] = fact;
			mpe[0][1] = mpe[1][0] = 0.0;
			mpe[2][2] = 2.0 * fact;
			mpe[0][2] = mpe[2][0] = mpe[1][2] = mpe[2][1] = 0.0;
		}
		// stiffness
		kpe[0][0] = kpe[1][1] = 0.0;
		kpe[0][1] = kpe[1][0] = 0.0;
		kpe[2][2] = 0.0;
		kpe[0][2] = kpe[2][0] = kpe[1][2] = kpe[2][1] = 0.0;
	}
}
