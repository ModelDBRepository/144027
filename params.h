#ifndef PARAMS_H
#define PARAMS_H

#include <math.h>
#include "lib.h"

class params
{
private:
	
public:


	// -------variables
	int	*display, ia_no, j_extra, neurons, newzap, norm, no_deps,
		no_display, no_mem, nu_no, seed, rezap, whozap;

	float	axons[2][2], *dneps, dnu, dt, dtzap, dw,
		dzap, eps, e_exc, e_inh, f, frac_inh,
		*ia_mean, *ia_sd, *ia_min, *ia_max, *ia_p,
		*nu_mean, *nu_sd, *nu_min, *nu_max, *nu_p, nu_zap[2],
		pspmax, pspmin, tau_m[2], tau_s[2], tmin, tmax,
		t_average, t_dzap_off, t_dzap_on, t_zap_off, t_zap_on, t_extra,
		vpsp[2][2], v_r, v_t, zap; 

	char	*suffix;

	// ---constructor
	params(int brgc, char** brgv, int prnt_cur, int prnt_in, int prnt_s);

	// ---get suffix
	void get_suffix(int brgc, char** brgv);

	// ---message explaining variables.
	void write_comline1();

	// ---sample list of variables.
	void write_comline2();

	// ---message explaining how current distribution is constructed.
	void write_comline3();

	// ---destructor
	~params()
	{
		delete [] display;
		delete [] ia_mean;
		delete [] ia_sd;
		delete [] ia_min;
		delete [] ia_max;
		delete [] ia_p;
	}
};

#endif
