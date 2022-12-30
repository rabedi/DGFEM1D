#include "globalTypesClasses.h"
#include "commonMacros.h"

void ElementProperties::Initialize_ElementProperties()
{
	c = sqrt(E / rho);
	Z = c * rho;
	time_e = hE / c;
}
