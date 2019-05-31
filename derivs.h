#ifndef DERIVS_H
#define DERIVS_H

#include "float2.h"
#include "network.h"

class derivs
{
private:
	int	neurons;
	float	*alpha, *deltai, *drho, *ihat, *mu, *rho, *tau;
	float	e_e, e_i, tau_s[2];

public:
	// number of variables
	int no_vars;

	// constructor
	derivs() {};

	// initialize
	void init_z(network& p, float2& z);

	// set parameters
	void set_params(network& p);

	// compute derivateves
	void func(float** dz, float2& z, float t);

	// destructor
	~derivs() {};

};

#endif
