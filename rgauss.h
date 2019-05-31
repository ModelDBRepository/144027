#ifndef GAUSS_H
#define GAUSS_H

//  headers
#include <cmath>
#include <cstdlib>
#include <lib.h>

//  definitions

//  prototypes
class rgauss
{
private:
	int n;
	float* y;

public:
	// constructor
	rgauss(int nn);

	// destructor
	~rgauss();

	// return random number with gaussian distribution.  second
	// call truncates from below, third trucates from below and above.
	float ran(float sigma);
	float ran(float sigma, float xmin);
	float ran(float sigma, float xmin, float xmax);

	// return array of numbers with gaussian distribution.  second
	// call truncates from below, third trucates from below and above.
	float* ran(float sigma, int ntot);
	float* ran(float sigma, int ntot, float xmin);
	float* ran(float sigma, int ntot, float xmin, float xmax);

	// return array of numbers with gaussian distribution.  second
	// call truncates from below, third trucates from below and above.
	// allow the possibility of saving alpha and beta for future calls.
	float* ran(
		float sigma, int ntot, float xmin, float xmax,
		float& alpha, float& beta);
};

#endif

