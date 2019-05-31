#include "network.h"

// setgauss needs to be included here because derivs has a #include network.h.
#include "setgauss.c"
#include <cstring>;

using namespace std;

/*
	Peter Latham
        May, 2002
	version 1

	---In a departure from the sn series, excitatory parameters
	   are first. So, for example, num[0] = number of excitatory
	   neurons and num[1] = number of inhibitory ones.

	---public variables
	deltav		- vthreshold-vrest, in mV.
	display		- array telling which neurons are displayed.
	dt		- time step in ms.
	dzap		- phase gets mulitplied by dezap, which must
			  be less than 1. this brings phase near zero
			  and turns neuron off.
	eta[i][j]	- eta[i] specifies which neurons are involved
			  in memory i. see set_mem for details.
	enorm[i]	- normalized reversal potential for cell of type i.
	f		- fraction of neurons involved in a memory.
	j_extra		- which extra spike appears.
	mu[i]		- normalized applied current for each cell.
	neurons		- number of neurons.
	no_display	- number of neurons that are displayed
	no_mem		- number of memories.
	num[i]		- number of neurons of type i.
	p_zap[i]	- probability of a zapping neuron firing per time step.
			  first arg is for on zapping; second is for off.
	suffix		- suffix that is appended to output files.
	tau[i]		- time constant of cell i, ms.
	tau_s[i]	- time constant for synapses of type i, ms.
	tmax		- maximum time of simulation, in ms.
	tmin		- starting time of simulation, in ms.
			  we write to file only for t >= 0, so setting
			  tmin negative will get rid of transients.
	type[i]		- type of neuron i. 0=exc, 1=inh.
	t_zap_on	- time in ms when zapping starts.
	t_zap_off       - time in ms when zapping starts.
	t_dzap_on	- time in ms when dezapping (hyperpolarizing) starts.
	t_dzap_off      - time in ms when dezapping (hyperpolarizing) stops.
	t_extra		- time at which extra spike appears.
	vbar		- (vrest+vthreshold)/2, in mV.
	w[j][i]		- weight, from neuron j to neuron wn[j][i].
	wn[j][i]        - index of nonzero weights.  e.g., if
			  wn[n][i] = j, then weight from neuron n to
			  neuron j is w[n][i].
	wntot[j]	- number of connections neuron j makes.
			  typical loop. z.x[0][n] is presynaptic,
			  z.x[2][wn[n][i]] is postsynaptic.

			  for (n=0; n < neurons; n++)
			  {
			     for (i=0; i < wntot[j]; i++)
			     {
			  	z.x[0][wn[n][i]] = w[n][i] ... z.x[0][n]
			     }
			  }
	who_zap[i]	- list of neurons within a sub-population.
			  who_zap[i]=1 means neuron is in sub-population,
			  who_zap[i]=0 means it is not.
	zap		- add phase zap to any neuron that is zapped.
*/

// ---constructor
network::network(int argc, char** argv, print pr)
{
	// ---get raw parameters
	params raw(argc, argv, pr.cur, pr.in, pr.s);


	// ---miscellaneous
	no_display = raw.no_display;
	display = (int*) calloc(no_display, sizeof(int));
	for (int i=0; i < no_display; i++) display[i] = raw.display[i];


	// ---parameters related to numbers of neurons.
	neurons	= raw.neurons;
	num[0]  = (int) (raw.frac_inh*neurons + 0.5) > neurons ? 0 :
			neurons - (int) (raw.frac_inh*neurons + 0.5);
	num[1]  = neurons-num[0];
	type    = set2(num[0], neurons, 0, 1);

	// ---firing rate quantization
	dnu     = raw.dnu;


	// ---parameters related to voltage.
	vbar	 = (raw.v_t + raw.v_r)/2.0;
	deltav	 = raw.v_t - raw.v_r;
	enorm[0] = (raw.e_exc-vbar)/deltav;
	enorm[1] = (raw.e_inh-vbar)/deltav;


	// ---time constants/interpoloation
	tau     = set2(num[0], neurons, raw.tau_m[0], raw.tau_m[1]);
	for (int i=0; i < 2; i++) tau_s[i] = raw.tau_s[i];


	// ---parameters related to time
	dt        = raw.dt;
	tmax      = raw.tmax;
	tmin      = raw.tmin;
	t_average = raw.t_average;


	// ---associative memory
	zap	   = raw.zap;
	dzap	   = raw.dzap;
	t_zap_on   = raw.t_zap_on;
	t_zap_off  = raw.t_zap_off;
	t_dzap_on  = raw.t_dzap_on;
	t_dzap_off = raw.t_dzap_off;
	no_mem     = raw.no_mem;
	f	   = raw.f;
	p_zap[0]   = raw.nu_zap[0]*(dt/1000);
	p_zap[1]   = raw.nu_zap[1]*(dt/1000);

	// ---rezap
	rezap      = raw.rezap;
	dtzap      = raw.dtzap;
	newzap     = raw.newzap;

	if (no_mem > 0)
	{
		deps = (float*) calloc(no_mem, sizeof(float));
		for (int i=0; i < no_mem; i++) deps[i]=1;
		for (int i=0; i < raw.no_deps/2; i++)
		{
			int m = (int) (raw.dneps[2*i] + 0.5) - 1;
			if (m < 0 || m >= no_mem)
				write_warning("n_deps=",m+1,"; out of range.");
			else deps[m] = raw.dneps[2*i+1];
		}
	}

	// ---extra spikes
	t_extra = raw.t_extra;
	j_extra = raw.j_extra;

	// ---output files
	suffix = (char*) calloc(1024, sizeof(char));
	strcpy(suffix, raw.suffix);


	// ---start random number generator
	srand(raw.seed);


	// ---gaussian random variable with 100*neurons possible values
	rgauss g(100*neurons);

	// ---applied current
	set_current(g, raw);

	// ---firing rate. quantized at dnu.
	nu_in = set_distribution(g, raw.nu_no, raw.nu_p, raw.nu_mean,
		raw.nu_sd, raw.nu_min, raw.nu_max, dnu);

	// ---weights. don't compute if pr.out=1.
	if (!pr.out)
	{
		float wmin, wmax;
		float** w_to_v = set_weights(g, raw, wmin, wmax);
	
		// ---memories
		if (no_mem)
		{
			// ---set memories
			set_mem(raw, wmin, wmax);

			// ---memories to zap
			who_zap = (int*) calloc(num[0], sizeof(int));
			if (raw.whozap >= 0) for (int i=0; i < num[0]; i++)
				who_zap[i] = eta[raw.whozap][i] > 0 ? 1 : 0;
			else for (int i=0; i < num[0]; i++)
				who_zap[i] = drand() < raw.f ? 1 : 0;
		}

		// ---print mean and standard deviation of weights. also
		//    compute actual max excitatory and inhibitory weights.
		float* vmax = get_stats(w_to_v);

		// ---diagnostic: print weights.
		if (pr.weights) write_weights(pr.weights, vmax, w_to_v);
	}

	// ---write parameters and exit if pr.out is on.
	if (pr.out) write_params(raw);
}

// ---set applied current, ia
void network::set_current(rgauss& g, params& raw)
{
/*
	The distribution of currents is the weighted sum of truncated
	Gaussian distributions, for mainly historical reasons.
	The mean, standard deviation, minimum and maximum values of the
	truncated Gaussian are ia_mean, ia_sd, ia_min, ia_max, and ia_sd,
	respectively. The number of weighted Gaussians is equal to the
	number of variables specified in ia_mean in the namelist.

	The variable that determines the weight of each distribution
	is ia_p. Unfortunately, this variable is a little bit complicated,
	because as an add-on I decided to let the excitatory and
	inhibitory neurons have a different distribution. So here's
	the convention:

		If all the ia_p are non-negative, then the excitatory and
		inhibitory neurons have the same distribution.

		If any of the ia_p are negative, then the excitatory neurons
		go with the distributions with non-negative ia_p and the
		inhibitory neurons go with the distributions with negative
		ia_p. For both types the total probability is forced to
		add to 1, and of course |ia_p| will be used to compute the
		probability.
	
	The code will exit if there are no distributions, or only one
	distribution with negative probability.
*/
	float* ia = set_distribution(g, raw.ia_no, raw.ia_p, raw.ia_mean,
		raw.ia_sd, raw.ia_min, raw.ia_max, 0.0);

	// ---now compute normalized drive.
	mu = (float*) calloc(neurons, sizeof(float));
	for (int n=0; n < neurons; n++) mu[n] = ia[n]/deltav - 0.25;

	// ---output excitatory and inhibitory current distribution
	FILE* ia_e = fopen(cat2("mu_e", raw.suffix), "w");
	for (int i=0; i < num[0]; i++)
		fprintf(ia_e, "%10.6f %10.6f\n", ia[i], mu[i]);
	fclose(ia_e);
	FILE* ia_i = fopen(cat2("mu_i", raw.suffix), "w");
	for (int i=num[0]; i < neurons; i++)
		fprintf(ia_i, "%10.6f %10.6f\n", ia[i], mu[i]);
	fclose(ia_i);

	free(ia);
}

// ---set distribution
float* network::set_distribution(rgauss& g, int n, float* p,
	float* mean, float* sd, float* min, float* max, float dx)
{
/*
	Constructs a distribution made up of a weighted sum of
	truncated Gaussian distributions, for mainly historical reasons.
	
	The code will exit if there are no distributions, or only one
	distribution with negative probability.

	INPUT:
	g	- object containing gaussian random numbers.
	n	- number of truncated Gaussian distributions. of 0, return
		  all 0's.
	p	- probability of each one occurring, with a twist:

		  If all the p are non-negative, then the excitatory
		  and inhibitory neurons have the same distribution.

		  If any of the p are negative, then the excitatory
		  neurons go with the distributions with non-negative
		  p and the inhibitory neurons go with the
		  distributions with negative ia_p. For both types the
		  total probability is forced to add to 1, and of
		  course |p| will be used to compute the probability.

		  For this to work, of course, this routine must have
		  access to the number of excitatory and inhibitory
		  neurons, which it does (through num[k]).

	mean	- mean of each distribution.
	sd	- standard deviation of each distribution.
	min	- min of each distribution.
	max	- max of each distribution.
	dx	- if > 0, values quantized by rounding to nearest integer*dx.

	RETURNS:
	x	- variable whose distribution is specified, as described
		  above. note that the excitatory and inhibitory portions
		  of x (the first num[0] and last num[1] elements) may have
		  different distributions.
*/
	// ---space for output
	float* x = (float*) calloc(neurons, sizeof(float));

	// ---if no data, return all zeros.
	if (n == 0) return x;

	// ---reserve (slightly too much) space
	float** prob = cfloat(2, n);
	float** f = cfloat(n, 5);

	// ---first, find probability distributions for each type.
	int i;
	for (i=0; i < n; i++) if (p[i] < 0) break;
	if (i == n)
	{
		// ---all distributions are non-negative.
		for (i=0; i < n; i++) for (int k=0; k < 2; k++)
			prob[k][i] = p[i];
	}
	else
	{
		// ---non-negative probabilities go to excitatory neurons,
		//    negative probabilities go to inhibitory neurons,
		for (i=0; i < n; i++) for (int k=0; k < 2; k++)
			if (p[i] < 0)
			{
				prob[1][i] = -p[i];
				prob[0][i] = 0;
			}
			else
			{
				prob[1][i] = 0;
				prob[0][i] = +p[i];
			}
	}

	/* ---f defined as follows:
	      f[i][j] - i labels distribution; it runs from 0 to n-1.
	                j=0: mean of gaussian.
	                  1: sd of gaussian.
	                  2: min of gaussian.
	                  3: max of gaussian.
                          4: probability that an element gets assigned to this
			     group.
	*/
	
	// ---set output
	float gmin=0;
	for (int k=0; k < 2; k++)
	{
		int cnt=0;
		for (int i=0; i < n; i++)
		{
			if (prob[k][i] > 0)
			{
				f[cnt][0] = mean[i];
				f[cnt][1] = sd[i];
				f[cnt][2] = min[i];
				f[cnt][3] = max[i];
				f[cnt][4] = prob[k][i];
				cnt++;

				// ---keep track of global min
				gmin = min[i] < gmin ? min[i] : gmin;
			}
		}
		float* x_tmp = setgauss(num[k], cnt, f, 0, g);

		int imin = k == 0 ? 0 : num[0];
		int imax = k == 0 ? num[0] : neurons;
		for (int i=imin; i < imax; i++) x[i] = x_tmp[i-imin];
		free(x_tmp);
	}

	// ---cleanup
	cfree(2, prob);
	cfree(n, f);

	// ---quantize if dx > 0
	if (dx > 0)
	{
		int imin = gmin >= 0 ? 0 : (int) (gmin/dx) - 1;
		for (int i=0; i < neurons; i++)
			x[i] = dx*(imin + (int) (x[i]/dx + 0.5 - imin));
	}

	// ---return
	return x;
}

// ---set weights, w
float** network::set_weights(rgauss& g, params& raw, float& wmin, float& wmax)
{
/*
	set weights.

	INPUT:
	g	- Gaussian random variable.
	raw	- object containing all sorts of parameters.

	OUTPUT:
	wmin	- minimum allowed weight; needed for clipping memories.
	wmax	- maximum allowed weight; needed for clipping memories.

	RETURNS:
	w_to_v	- factor used to translate from weights to psp size. formula is:

		     vpsp[type[i][type[j]] =
		        w_to_v[type[i]][type[wn[i][j]]]*w[i][j]

		  remember, convention is w[i][j] is weight from neuron i
		  to neuron wn[i][j]. w_to_v[i][j] is from neuron of type i to
		  neuron of type j.

	COMPUTES:
	w[i][j]	- weight, from neuron i to neuron wn[i][j].
	wn[i][j]- index of postsynaptic neuron.
	wntot[i]- number of connections neuron i makes.
*/

	// ---probability of connection. pconnect[i][j] is the probability
	//    that a neuron of type j connects to a neuron of type i.
	float** pconnect = cfloat(2, 2);
        for (int i=0; i < 2; i++) for (int j=0; j < 2; j++)
                pconnect[i][j] = raw.axons[j][i]/neurons;

	// ---normalized PSP sizes, from neuron i to neuron j.
	float** v_to_w = cfloat(2, 2);
	float** w_to_v = cfloat(2, 2);
	float tau_m[2]    = {raw.tau_m[0], raw.tau_m[1]};
	for (int i=0; i < 2; i++) for (int j=0; j < 2; j++)
	{
		float dum = tau_m[j] == tau_s[i] ? 1 :
			log(tau_m[j]/tau_s[i])/(tau_m[j]/tau_s[i]-1);

		// ---vpsp[i][j] is from neuron i to neuron j.
		v_to_w[i][j] = (raw.vpsp[i][j]/deltav)
			  * (1/(enorm[i]+0.5))
			  * (tau_m[j]/tau_s[i])*exp(dum);
		w_to_v[i][j] = raw.vpsp[i][j]/v_to_w[i][j];

		// ---w_to_v is positive by convention
		w_to_v[i][j] *= w_to_v[i][j] > 0 ? 1 : -1;

		// ---set wbar, positive by convention
		wbar[i][j] = v_to_w[i][j] > 0 ? v_to_w[i][j] : -v_to_w[i][j];
	}

	// ---minimum and maximum weights between excitatory neurons.
	//    needed for setting memories.
	wmin = v_to_w[0][0]*raw.pspmin/raw.vpsp[0][0];
	wmax = v_to_w[0][0]*raw.pspmax/raw.vpsp[0][0];

	// ---set weights. w[j][i] is weight from neuron j to neuron wn[j][i]
	w = new float*[neurons];
	wn = new int*[neurons];
	wntot = new int[neurons];
	float sqrt3 = sqrt(3);
	for (int j=0; j < neurons; j++)
	{
		ostrstream w_ptr;
		ostrstream wn_ptr;
		wntot[j]=0;
		for (int i=0; i < neurons; i++)
		{
			if (drand() < pconnect[type[i]][type[j]] && i != j)
			{
				float ww = v_to_w[type[j]][type[i]]*
					(1+raw.dw*drand(-sqrt3, sqrt3));
				w_ptr.write((char*) &ww, sizeof(ww));
				wn_ptr.write((char*) &i, sizeof(i));
				wntot[j]++;
			}
		}

		// freeze ostrstream
		w[j] = (float*) w_ptr.str();
		wn[j] = (int*) wn_ptr.str();
	}

	return w_to_v;
}

// ---set memories
void network::set_mem(params& raw, float wmin, float wmax)
{
/*
	Set memories. For each memory, the connection strength from
	neuron j to neuron i is increased by

		eps * (eta[i] + f*norm) * eta[j]
	
	where only excitatory neurons are involved and:

		eps	= strength of memory normalized to sqrt(p*f*(1-f))
			  where p is the number of memories and f is
			  the fraction of neurons involved in a
			  memory (see eta).
		eta	= 1-f with probability f and -f with probability 1-f.
		norm	= 0: both pre and post-synaptic normalization,
			  1: postsynatpic normalization only.

	After modifying the connection strength, the weights are
	clipped: they can't fall below wmin or go above wmax.
*/
	// ---space for eta.
	if (no_mem) eta = cfloat(no_mem, num[0]);

	// ---modify memories
	for (int m=0; m < no_mem; m++)
	{
		// ---normalize memory strength.
		float wee = deps[m]*raw.eps/sqrt(raw.f*(1-raw.f));

		// --create eta.
		set_eta(m, raw.f);

		// ---modify connection strengths among excitatory neurons,
		//    from neuron j to neurons wn[j][i].
		for (int j=0; j < num[0]; j++)
		{
			for (int i=0; i < wntot[j]; i++) if (!type[wn[j][i]])
			{
				w[j][i] += wee*(eta[m][wn[j][i]]+
				        raw.f*raw.norm)*eta[m][j];
			}
		}
	}

	// ---clip memories.
	for (int j=0; j < neurons; j++) if (type[j] == 0)
	{
		for (int i=0; i < wntot[j]; i++)
		{
			w[j][i] = w[j][i] < wmin ? wmin : w[j][i];
			w[j][i] = w[j][i] > wmax ? wmax : w[j][i];
		}
	}

	// ---make list of which memories each neuron participates in.
	//    wmem stands for "which memory".
	int* wtmp = (int*) calloc(no_mem, sizeof(int));
	nwmem = (int*) calloc(num[0], sizeof(int));
	wmem = (int**) calloc(num[0], sizeof(int*));
	for (int i=0; i < num[0]; i++)
	{
		for (int j=0; j < no_mem; j++) if (eta[j][i] > 0)
			wtmp[nwmem[i]++]=j;
		wmem[i] = (int*) calloc(nwmem[i], sizeof(int));
		for (int j=0; j < nwmem[i]; j++) wmem[i][j] = wtmp[j];
	}

	free(wtmp);
}

// ---create memory vector, eta.
void network::set_eta(int m, float f)
{
/*
	set eta:  if neuron is excitatory:  1-f with probability f
					    -f  with probability 1-f
		  if neuron is inhibitory:  no value
*/
	for (int i=0; i < num[0]; i++)
		eta[m][i] = (drand() < f && !type[i]) ? 1-f : -f;
}

// ---get statistics on memory:  mean and variance
float* network::get_stats(float** w_to_v)
{
	// ---reserve space
	double** mean = init_constant(2, 2, (double) 0);
	double** sigma = init_constant(2, 2, (double) 0);
	int** cnt = init_constant(2, 2, 0);

	// ---means, min, max
	float* vmax = (float*) calloc(2, sizeof(float));
	for (int i=0; i < neurons; i++) for (int j=0; j < wntot[i]; j++)
	{
		// ---index of the postsynaptic neuron
		int k = wn[i][j];

		// ---first and second moments
		mean[type[i]][type[k]] += w[i][j];
		sigma[type[i]][type[k]] += w[i][j]*w[i][j];
		cnt[type[i]][type[k]]++;

		// ---min and max
		double mv = w_to_v[type[i]][type[wn[i][j]]]*w[i][j];
		vmax[type[i]] = mv > vmax[type[i]] ? mv : vmax[type[i]];
	}

	// ---normalize and compute standard deviation
	for (int k=0; k < 2; k++) for (int l=0; l < 2; l++)
	{
		mean[k][l] /= cnt[k][l] == 0 ? 1 : cnt[k][l];
		sigma[k][l] /= cnt[k][l] == 0 ? 1 : cnt[k][l];
		sigma[k][l] = sigma[k][l] - mean[k][l]*mean[k][l];
		sigma[k][l] = sigma[k][l] < 0 ? 0 : sqrt(sigma[k][l]);
	}

	// ---print results
	cerr << "           mean       standard deviation" << endl;
	fprintf(stderr, "E --> E: %10.6f %10.6f\n", mean[0][0], sigma[0][0]);
	fprintf(stderr, "E --> I: %10.6f %10.6f\n", mean[0][1], sigma[0][1]);
	fprintf(stderr, "I --> E: %10.6f %10.6f\n", mean[1][0], sigma[1][0]);
	fprintf(stderr, "I --> I: %10.6f %10.6f\n", mean[1][1], sigma[1][1]);

	return vmax;
}

// ---diagnostic: print weights
void network::write_weights(int nbins, float* vmax, float** w_to_v)
{
	// ---output weights for display neurons.
	for (int i=0; i < no_display; i++)
	{
		// ---presynaptic neuron
		int n = display[i];

		// ---open file
		char s[1024];
		sprintf(s, "w_%d%s", n, suffix);
		FILE* w_out = fopen(s, "w");

		// ---write
		for (int j=0; j < wntot[n]; j++)
			fprintf(w_out, "%d %10.6f\n", wn[n][j],
				w_to_v[type[n]][type[wn[n][j]]]*w[n][j]);

		// ---close file
		fclose(w_out);
	}

	// ---make histograms. first step is range of weights. weights
	//    from inhibitory and excitatory cells have different
	//    signs, so they're binned differently.

	// ---reserve space
	float** wbins = cfloat(2, nbins);
	int** counts = cint(3, nbins);

	// ---expand max weights so that we don't have spillover.
	for (int i=0; i < 2; i++) vmax[i] *= 1.0001;

	// ---delta w
	float delta[2] = {vmax[0]/nbins, vmax[1]/nbins};

	// ---weights for the bins.
	for (int k=0; k < 2; k++)
		for (int i=0; i < nbins; i++) wbins[k][i] = (i+0.5)*delta[k];

	// ---counts
	for (int i=0; i < neurons; i++)
	{
		// ---determine k:
		//       0 for excitatory non-memory, 
		//       1 for excitatory memory.
		//       2 for inhibitory,
		int k = type[i] == 1 ? 2 :
			(no_mem == 0 ? 0 : (eta[0][i] < 0 ? 0 : 1));

		// ---get counts
		for (int j=0; j < wntot[i]; j++)
			counts[k][(int)
				(w_to_v[type[i]][type[wn[i][j]]]*
				 	w[i][j]/delta[type[i]])]++;
	}

	// ---space for filename
	char s[1024];

	// ---output excitatory histogram
	sprintf(s, "w_e%s", suffix);
	FILE* w_exc = fopen(s, "w");
	for (int i=0; i < nbins; i++)
		fprintf(w_exc, "%10.6f %d %d\n",
			wbins[0][i], counts[0][i], counts[1][i]);
	fclose(w_exc);

	// ---output inhibitory histogram
	sprintf(s, "w_i%s", suffix);
	FILE* w_inh = fopen(s, "w");
	for (int i=0; i < nbins; i++)
		fprintf(w_inh, "%10.6f %d\n", wbins[1][i], counts[2][i]);
	fclose(w_inh);

	// ---free memory
	cfree(2, wbins);
	cfree(3, counts);
}

// ---set values
int* network::set2(int n1, int n2, int p1, int p2)
{
/*
	set first n1 values to p1, next n2-n1 values to p2.
*/
	int* x = (int*) calloc(n2, sizeof(int));
	for (int i=0 ; i < n1; i++) x[i]=p1;
	for (int i=n1; i < n2; i++) x[i]=p2;
	return x;
}

// ---set values
float* network::set2(int n1, int n2, float p1, float p2)
{
/*
	set first n1 values to p1, next n2-n1 values to p2.
*/
	float* x = (float*) calloc(n2, sizeof(float));
	for (int i=0 ; i < n1; i++) x[i]=p1;
	for (int i=n1; i < n2; i++) x[i]=p2;
	return x;
}

// ---write parameters
void network::write_params(params& raw)
{
	cout << "neurons:             " << neurons       << endl;
	cout << "exc. neurons:        " << num[0]        << endl;
	cout << "inh. neurons:        " << num[1]        << endl;
	cout << "dt:                  " << dt            << endl;
	cout << "tmin:                " << tmin          << endl;
	cout << "tmax:                " << tmax          << endl;
	cout << "exc. tau_s:          " << tau_s[0]      << endl;
	cout << "inh. tau_s:          " << tau_s[1]      << endl;
	cout << "exc. norm v_reverse  " << enorm[0]      << endl;
	cout << "inh. norm v_reverse  " << enorm[1]      << endl;

	// ---PSP sizes
	for (int i=0; i < 2; i++) for (int j=0; j < 2; j++) cout
		<< "PSP  ("
		<< (i ? "I" : "E") << " --> "
		<< (j ? "I" : "E") << ")       "
		<< " "
		<< raw.vpsp[i][j]
		<< endl;

	cout << "no_display:         " << no_display << endl;
	cout << "displayed neurons:  ";
	for (int i=0; i < no_display; i++) cout << display[i] << " ";
	cout << endl;

	char s[1024];
	if (!strcmp(suffix, "")) strcpy(s, "");
	else strcpy(s, &suffix[1]);
	cout << "suffix:             " << s          << endl;

	if (no_mem)
	{
		cout << endl;
		cout << "no_mem:             " << no_mem     << endl;
		cout << "zap:                " << zap        << endl;
		cout << "dzap:               " << dzap       << endl;
		cout << "f:                  " << raw.f      << endl;
		cout << "norm:               " << raw.norm   << endl;
		cout << "whozap:             " << raw.whozap << endl;
		cout << "nu_zap:             " << raw.nu_zap[0] << " "
					       << raw.nu_zap[1] << endl;
		cout << "t_zap_on:           " << t_zap_on   << endl;
		cout << "t_zap_off:          " << t_zap_off  << endl;
		cout << "t_dzap_on:          " << t_dzap_on  << endl;
		cout << "t_dzap_off:         " << t_dzap_off << endl;

		if (raw.no_deps > 0)
		{
			cout << "modified memories:";
			int notfirst=0;
			for (int i=0; i < raw.no_deps/2; i++)
			{
				int m = (int) (raw.dneps[2*i] + 0.5) - 1;
				if (m >= 0 && m < no_mem)
				{
					if (notfirst++)
						cout << "                  ";
					fprintf(stdout, "%3d  %5.3f\n",
						m, raw.dneps[2*i+1]);
				}
			}
		}
	}

	exit(1);
}
