/*
	cumulative probability for poisson distribution. returns
	a variable q[i] such that prob(q[i]=k) = p^k exp(-p)/k!.
	q[i] has n elements, ordered from smallest to largest.
*/

// headers
#include <math.h>
#include <stdlib.h>

double* dpoisson(float p, int n)
{
/*
	INPUT:
	p	- probability.
	n	- number of points in q: q runs from 0 to n-1.

	RETURNS:
	q[i]	- probability that q[i]=k is p^k exp(-p)/k!. q[i] is ordered
		  from smallest to largest.
*/
	// ---reserve space
	double* q = (double*) calloc(n, sizeof(double));

	// ---set everybody to zero if p=0.
	if (p <= 0)
	{
		for (int i=0; i < n; i++) q[i]=0.0;
		return q;
	}

	// ---constants
	double logp = log(p);

	// --- up
	int j1, j2=0;
	double  logfac=0, qacc=0;
	for (int i=0;; i++)
	{
		double q0 = exp(i*logp - p - logfac);

		j1=j2;
		j2=(int) (n*(qacc+q0));
		for (int j=j1; j < j2; j++) q[j]=i;

		qacc += q0;
		logfac += log(1+i);

		if (j2 == n-1) break;
	}
	q[n-1]=q[n-2]+1;

	return q;
}
