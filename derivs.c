#include "derivs.h"

/*
	Peter Latham
	May, 2002
	version 1
*/

// ---initialize
void derivs::init_z(network& p, float2& z)
{
/*
	initialize z.  if variables are added, be sure to change no_vars.

	INPUT is class p, from which following members are used:
	enorm[i]- normalized reversal potential for cells of type i.
	mu[i]	- normalized applied current.
	neurons	- number of neurons.
	tau[i]	- cell time constant, ms.
	tau_s[i]- time constant for synapses of type i, ms.

	OUTPUT:
	z		- dynamical variable, plus information about itself.

	SETS:
	no_vars		- number of variables.
*/
	// ---change this whenever the number of variables changes.
	no_vars = 3;

	// ---set z.n
	z.n = no_vars;

	// ---set z.m
	z.m = new int[no_vars];
	for (int n=0; n < no_vars; n++) z.m[n] = p.neurons;

	// ---set variables. phase is initially random.
	z.x = newdouble(no_vars, p.neurons);
	float pi = 4.0*atan(1.0);
	for (int n=0; n < p.neurons; n++)
	{
		z.x[0][n] = drand(-pi, 0);	// phase
		z.x[1][n] = 0;			// excitatory drive (I_E)
		z.x[2][n] = 0;			// inhibitory drive (I_I)
	}

	// ccc use z.x[0][n] = -pi*fmod(n_ic*drand(), 1.0)
}

// ---set parameters
void derivs::set_params(network& p)
{
/*
	INPUT:
	p	- class containing single neuron properties.

*/

	// ---reserve space
	mu     = (float*) calloc(p.neurons, sizeof(float));
	tau    = (float*) calloc(p.neurons, sizeof(float));

	// ---set parameters
	e_e	= p.enorm[0];
	e_i	= p.enorm[1];
	neurons = p.neurons;

	for (int i=0; i < p.neurons; i++)
	{
		mu[i]     = p.mu[i];
		tau[i]    = p.tau[i];
	}
	for (int i=0; i < 2; i++) tau_s[i] = p.tau_s[i];
}

// ---compute derivatives
void derivs::func(float** dz, float2& z, float t)
{
/*
	compute derivatives needed by integrator

	INPUT:
	z.x[0]	- phase
	z.x[1]	- excitatory drive (I_E)
	z.x[2]	- inhibitory drive (I_I)
	t	- time

	OUTPUT:
	dz	- time derivatives

	PRIVATE VARIABLES:
	e_e	- normalized excitatory reversal potential.
	e_i	- normalized inhibitory reversal potential.
	mu[i]	- normalized applied current.
	neurons	- number of neurons.
	tau[i]	- cell time constant, ms.
	tau_s[i]- time constant for synapses of type i, ms.
*/
	// ---save cos and sin
	static double* cs = (double*) calloc(neurons, sizeof(double));
	static double* sn = (double*) calloc(neurons, sizeof(double));
	for (int i=0; i < neurons; i++)
	{
		cs[i] = cos(z.x[0][i]);
		sn[i] = sin(z.x[0][i]);
	}

	// ---phase. two possible nonlinearities
	for (int i=0; i < neurons; i++) dz[0][i] =
		(1 - cs[i] + mu[i]*(1+cs[i])
		 - (sn[i] - e_e*(1+cs[i]))*z.x[1][i]
		 - (sn[i] - e_i*(1+cs[i]))*z.x[2][i])/tau[i];

	// ---excitatory and inhibitory currents
	for (int i=0; i < neurons; i++) dz[1][i] = -z.x[1][i]/tau_s[0];
	for (int i=0; i < neurons; i++) dz[2][i] = -z.x[2][i]/tau_s[1];
}
