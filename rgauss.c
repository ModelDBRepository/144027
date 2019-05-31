/*


	Peter Latham
	July, 1997

*/

#include "rgauss.h"

// constructor
rgauss::rgauss(int nn)
{
	// set variables and reserve space
	this->n = nn;
	this->y = new float[n];

	/* divide gaussian into n intervals of equal probability.
	estimate area by expanding gaussian through quadratic.
	width in y of section with probability r is given approximately by

	dy = (r/f) - (f'/2f) * (r/f)^2 +[(f'/f)^2/2 - f''/6f)] * (r/f)^3

	where   dy = y-y0, f=f(y0), f'=df(y0)/dy0, f''=d^f(y0)/dy0^2,
		f is gaussian probability distribution. */

	int n1, n2;
	// sqrt2pi = sqrt(2*pi)
	float r, sqrt2pi = 2.5066283, y0;
	r=1/((float) n);

	// determine where to start.  that depends on whether n is odd or even.
	if (n%2 == 0)
	{
		y0 = 0;
		n1=n/2;
		n2=n/2-1;
	}
	else
	{
		float dum = 0.5*r*sqrt2pi;
		y0 = dum*(1+dum*dum/6);
		y[n/2] = 0;
		n1=n/2+1;
		n2=n/2-1;
	}

	for (int i=0; i < n/2; i++)
	{
		float dum = r*sqrt2pi*exp(y0*y0/2);
		float dy = dum*(1+dum*(0.5*y0+dum*(0.5+y0*y0)/3));
		y[n1+i] = y0+0.5*dy;
		y[n2-i] = -(y0+0.5*dy);
		y0 += dy;
	}
}

// destructor
rgauss::~rgauss()
{
	delete [] y;
}


// overloaded function.
float rgauss::ran(float sigma)
{
/*
	return random number with gaussian distribution.

	INPUT:
	sigma	- standard deviation of Gaussian.

	RETURNS:
	gaussian distributed random number with standard deviation sigma.
*/
	return sigma*y[(int) (n*drand())];
}

float rgauss::ran(float sigma, float xmin)
{
/*
	return random number with gaussian distribution truncated below.

	INPUT:
	sigma	- standard deviation of Gaussian.
	xmin	- minimum cutoff of truncated gaussian.

	RETURNS:
	truncated (from below) gaussian distributed random number
	with standard deviation (before truncation) sigma.
*/
	int i;

	// find cutoff
	for (i=0; i < n; i++) if (sigma*y[i] >= xmin) break;
	if (i >= n) i=n-1;
	float alpha = ((float) i) / ((float) n);

	float x = sigma*y[(int) (n*(alpha + (1-alpha)*drand()))];

	return x;
}

float rgauss::ran(float sigma, float xmin, float xmax)
{
/*
	return random number with gaussian distribution truncated below
	and above.

	INPUT:
	sigma	- standard deviation of Gaussian.
	xmin	- minimum cutoff of truncated gaussian.
	xmax	- maximum cutoff of truncated gaussian.

	RETURNS:
	truncated (from below and above) gaussian distributed random
	number with standard deviation (before truncation) sigma.
*/
	// check
	if (xmax < xmin)
	{
		float tmp = xmin;
		xmin = xmax;
		xmax = tmp;
	}
	int i, j;

	// find minimum cutoff
	for (i=0; i < n; i++) if (sigma*y[i] >= xmin) break;
	if (i >= n) i=n-1;

	// find maximum cutoff
	for (j=n-1; j >= 0; j--) if(sigma*y[j] <= xmax) break;
	if (j < 0) j=0;

	// check
	if (i == j)
	{
		if (i == 0) j=1;
		else i = j-1;
	}

	float alpha = ((float) i) / ((float) n);
	float beta = ((float) j) / ((float) n);

	float x = sigma*y[(int) (n*(alpha + (beta-alpha)*drand()))];

	return x;
}

float* rgauss::ran(float sigma, int ntot)
{
/*
	return array of random numbers with gaussian distribution.

	INPUT:
	sigma	- standard deviation of Gaussian.
	ntot	- total number of random numbers to return.

	RETURNS:
	array of gaussian distributed random numbers with standard
	deviation sigma.
*/
	float* x = new float[ntot];
	for (int i=0; i < ntot; i++) x[i] = sigma*y[(int) (n*drand())];
	return x;
}

float* rgauss::ran(float sigma, int ntot, float xmin)
{
/*
	return array of random numbers with gaussian distribution truncated
	below.

	INPUT:
	sigma	- standard deviation of Gaussian.
	ntot	- total number of random numbers to return.
	xmin	- minimum cutoff of truncated gaussian.

	RETURNS:
	array of truncated (from below) gaussian distributed random numbers
	with standard deviation (before truncation) sigma.
*/
	int i;

	// find cutoff
	for (i=0; i < n; i++) if (sigma*y[i] >= xmin) break;
	if (i >= n) i=n-1;
	float alpha = ((float) i) / ((float) n);

	float* x = new float[ntot];
	for (i=0; i < ntot; i++)
		x[i] = sigma*y[(int) (n*(alpha + (1-alpha)*drand()))];

	return x;
}

float* rgauss::ran(float sigma, int ntot, float xmin, float xmax)
{
/*
	return array of random numbers with gaussian distribution truncated
	below and above.

	INPUT:
	sigma	- standard deviation of Gaussian.
	ntot	- total number of random numbers to return.
	xmin	- minimum cutoff of truncated gaussian.
	xmax	- maximum cutoff of truncated gaussian.

	RETURNS:
	array of truncated (from below and above) gaussian distributed random
	numbers with standard deviation (before truncation) sigma.
*/

	// check
	if (xmax < xmin)
	{
		float tmp = xmin;
		xmin = xmax;
		xmax = tmp;
	}
	int i, j;
	
	// find minimum cutoff
	for (i=0; i < n; i++) if (sigma*y[i] >= xmin) break;
	if (i >= n) i=n-1;

	// find maximum cutoff
	for (j=n-1; j >= 0; j--) if(sigma*y[j] <= xmax) break;
	if (j < 0) j=0;
	
	// check
	if (i == j)
	{
		if (i == 0) j=1;
		else i = j-1;
	}

	float alpha = ((float) i) / ((float) n);
	float beta = ((float) j) / ((float) n);

	float* x = new float[ntot];
	for (i=0; i < ntot; i++)
		x[i] = sigma*y[(int) (n*(alpha + (beta-alpha)*drand()))];
	return x;
}

float* rgauss::ran(
	float sigma, int ntot, float xmin, float xmax,
	float& alpha, float& beta)
{
/*
	return array of random numbers with gaussian distribution truncated
	below and above.

	INPUT:
	sigma	- standard deviation of Gaussian.
	ntot	- total number of random numbers to return.
	xmin	- minimum cutoff of truncated gaussian.
	xmax	- maximum cutoff of truncated gaussian.
	alpha, beta-
		- if both these variables are nonnegative, then they are
		  used as cutoff values instead of xmin and xmax.
		  they are used for repeated calls to gauss with
		  the same xmin and xmax, since it can take a long
		  time to compute alpha and beta.

	OUTPUT:
	alpha, beta
		- returned for use on next call, if desired.  see above.

	RETURNS:
	array of truncated (from below and above) gaussian distributed random
	numbers with standard deviation (before truncation) sigma.
*/
	float* x = new float[ntot];

	if (alpha < 0 || beta < 0)
	{
		// ---compute alpha and beta from xmin and xmax

		// check
		if (xmax < xmin)
		{
			float tmp = xmin;
			xmin = xmax;
			xmax = tmp;
		}
		int i, j;
	
		// find minimum cutoff
		for (i=0; i < n; i++) if (sigma*y[i] >= xmin) break;
		if (i >= n) i=n-1;

		// find maximum cutoff
		for (j=n-1; j >= 0; j--) if(sigma*y[j] <= xmax) break;
		if (j < 0) j=0;
	
		// check
		if (i == j)
		{
			if (i == 0) j=1;
			else i = j-1;
		}

		alpha = ((float) i) / ((float) n);
		beta = ((float) j) / ((float) n);
	}

	// generate random numbers
	for (int i=0; i < ntot; i++)
		x[i] = sigma*y[(int) (n*(alpha + (beta-alpha)*drand()))];

	return x;
}
