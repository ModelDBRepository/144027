#ifndef NETWORK_H
#define NETWORK_H

#include <math.h>
#include "rgauss.h"
#include "params.h"
#include "lib.h"

// ---structures
struct print
{
	int cur, f, in, out, s, weights;

};

class network
{
public:
	int	*display, j_extra, neurons, newzap, no_display, no_mem, num[2],
		*nwmem, rezap, *type, **wmem, *who_zap, **wn, *wntot;

	float	deltav, *deps, dnu, dzap, dt, dtzap,
		**eta, enorm[2], f, *mu, *nu_in, p_zap[2],
		*tau, tau_s[2], tmax, tmin, t_average, t_dzap_off,
		t_dzap_on, t_zap_off, t_zap_on, t_extra,
		vbar, **w, wbar[2][2], zap; 

	char	*suffix;

	// ---constructor
	network(int argc, char** argv, print pr);

	// ---destructor
	~network()
	{
		free(type);
		free(mu);
		free(tau);

		for (int i=0; i < neurons; i++)
		{
			delete [] w[i];
			delete [] wn[i];
		}
		delete [] display;
		delete [] w;
		delete [] wn;
		delete [] wntot;

		if (no_mem)
		{
			free(who_zap);
			cfree(no_mem, eta);
			free(nwmem);
			cfree(num[0], wmem);
		}
	}

private:
	// ---set applied current, ia
	void set_current(rgauss& g, params& raw);

	// ---set distribution
	float* set_distribution(rgauss& g, int n, float* p,
		float* mean, float* sd, float* min, float* max, float dx);

	// ---set weights, w
	float** set_weights(rgauss& g, params& raw, float& wmin, float& wmax);

	// ---set memories
	void set_mem(params& raw, float wmin, float wmax);

	// ---create memory vector, eta.
	void set_eta(int m, float f);

	// ---get statistics on memory:  mean and variance
	float* get_stats(float** w_to_v);

	// ---diagnostic: print weights
	void write_weights(int nbins, float* maxw, float** w_to_v);

	// ---set paraemters to one of two values.
	int* set2(int n1, int n2, int p1, int p2);
	float* set2(int n1, int n2, float p1, float p2);

	// ---write parameters
	void write_params(params& raw);
};

#endif
