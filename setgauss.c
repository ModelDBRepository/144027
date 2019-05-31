/*
	returns an array of variables with probability distribution

		p(x) = sum_j rho_j exp(-(mu_j-x_j)^2/2*sigma_j^2)/Z
	
	where the Gaussian may be truncated and Z is the normalization,
	chosen such that Z=int dx_xa^xb exp(-(mu_j-x_j)^2/2*sigma_j^2) where
	xa and xb are the truncation values. the rho_j sum to 1.

	1. make this float* rather than double*.
	2. add seed, but don't call srand unless seed > 0.
	3. allow nsplit=0.

*/

// headers
#include <math.h>
#include <stdlib.h>

float* setgauss(int n, int nf, float** f, int seed, rgauss& g)
{
/*
	n elements are assigned randomly to nf groups.

	f[i][4] gets normalized. if any numbers are negative, we just take
	the absolute value.

	INPUT:
	n	- total number of elements. output is size n.
	nf	- number of probability distributions.
	f[i][j]	- i labels distribution; it runs from 0 to nf-1.
		  j=0: mean of gaussian.
		    1: sd of gaussian.
		    2: min of gaussian.
		    3: max of gaussian.
		    4: probability that an element gets assigned to this group.
	seed	- if > 0, call srand(seed)

	RETURNS:
	x[i]	- elements with (sum of gaussian) distribution, as
		  discussed above.
*/

	// ---start random number generator if seed > 0.
	if (seed > 0) srand(seed);

	// ---muck with f
	double* facc = (double*) calloc(nf, sizeof(double));
	for (int i=0; i < nf; i++)
	{
		// ---normalize f[i][4]
		double pacc=0;
		for (int i=0; i < nf; i++)
		{
			f[i][4] = f[i][4] < 0 ? -f[i][4] : f[i][4];
			pacc += f[i][4];
		}
		for (int i=0; i < nf; i++) f[i][4] /= pacc;

		// ---make cumulative distribution
		facc[nf-1]=1.0;
		for (int i=nf-1; i > 0; i--)
			facc[i-1] = facc[i] - f[i][4];
	}

	// ---find out which element is in which group.
	int* dist = (int*) calloc(n, sizeof(int));
	int* cnt = (int*) calloc(nf, sizeof(int));
	for (int i=0; i < n; i++)
	{
		// ---assign neuron to distribution
		double rnum = drand();
		int j;
		for (j=0; j < nf; j++) if (rnum <= facc[j]) break;

		dist[i]=j;
		cnt[j]++;
	}

	// ---get random numbers
	float** xtmp = (float**) calloc(nf, sizeof(float*));
	for (int i=0; i < nf; i++) xtmp[i] =
		g.ran(f[i][1], cnt[i], f[i][2]-f[i][0], f[i][3]-f[i][0]);

	// ---set x
	float* x = (float*) calloc(n, sizeof(float));
	for (int i=0; i < nf; i++) cnt[i]=0;
	for (int i=0; i < n; i++)
		x[i] = f[dist[i]][0]+xtmp[dist[i]][cnt[dist[i]]++];

	free(cnt);
	free(dist);
	free(facc);
	cfree(nf, xtmp);
	return x;
}
