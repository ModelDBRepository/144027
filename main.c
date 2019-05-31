/*

	Peter Latham
	September, 2005
	version 1

	type "theta1d -h" or "theta1d" for details.

	ttd:
	- get rid of zap (if it really doesn't do anything).

	- make file containing rate at which +- pi is crossed illegally.
	- give whozap a random component.

	changes:
	09/05	- killed interp and nonlinearity.
		  removed bias (old conventions).

		  added t_extra and j_extra -- time of extra spike and
		  which neuron it occurs on, respectively.

*/

// ---headers
#include <cstdlib>
#include <fstream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cctype>
#include <iostream>

#include "derivs.h"
#include "network.h"
#include "params.h"
#include "float2.h"
#include "runge_kutta.h"
#include "lib.h"
#include "dpoisson.c"

// ---parameters

// ---prototypes
float* integrate(network& p, int write_on);
void spike(network& p, float t, float2& z, int& n_spiked, int* nspikes,
	int* who_spiked, int* cnterr);
void update_post(network& p, float2& z, int* nspikes, int* cnterr, int j);
double** input_dist(network& p, int* l, int& n, int n_dp);
void write(int* display, int n_spiked, int no_display, float t,
	int* who_spiked, int* count_total, int* aspikes, float2& z, network& p);
void count_ea(network& p);
void write_comline();
void write_filelist();

int main(int argc, char** argv)
{
	// ---command line stuff
	if (argc == 1) write_comline();
	if (argc > 1) if (!strcmp(argv[1], "-h")) write_comline();
	if (parse_comline(argc, argv, "-h", 0, 0)) write_comline();

	// ---get command line flags.
	struct print pr;
	pr.cur       = parse_comline(argc, argv, "-c", 0, 0);
	pr.f         = parse_comline(argc, argv, "-f", 0, 0);
	pr.in        = parse_comline(argc, argv, "-i", 0, 0);
	pr.out       = parse_comline(argc, argv, "-o", 0, 0);
	pr.s         = parse_comline(argc, argv, "-s", 0, 0);
	pr.weights   = parse_comline(argc, argv, "-w", 1, 0);
	int write_on = 1-parse_comline(argc, argv, "-nowrite", 0, 0);


	// ---give list of files if -f flag is present
	if (pr.f) write_filelist();

	// ---generate parameters. see network.c for details.
	network p(argc, argv, pr);

	// ---count number of endogenously active cells and print results
	count_ea(p);

	// ---output memories, if there are any
	if (p.no_mem)
	{
		ofstream out_mem(cat2("eta", p.suffix));
		for (int j=0; j < p.num[0]; j++)
		{
			for (int m=0; m < p.no_mem; m++)
				out_mem << p.eta[m][j] << " ";
			out_mem << endl;
		}
	}

	// ---integrate equations of motion; return average firing rate
	float* nu = integrate(p, write_on);

	// write to file fc
	ofstream out_fc_e(cat2("fc_e", p.suffix));
	ofstream out_fc_i(cat2("fc_i", p.suffix));
	for (int i=0; i < p.num[0]; i++) out_fc_e << nu[i] << endl;
	for (int i=p.num[0]; i < p.neurons; i++) out_fc_i << nu[i] << endl;

	// ---write average firing inhibitory and excitatory rates to file f
	float nubar[2];  nubar[0] = nubar[1] = 0;
	for (int i=0; i < p.neurons; i++) nubar[p.type[i]] += nu[i];
	for (int i=0; i < 2; i++) nubar[i] /= (p.num[i] > 0 ? p.num[i] : 1);
	ofstream out_f(cat2("f", p.suffix));
	out_f << nubar[0] << " " << nubar[1] << endl;
}

float* integrate(network& p, int write_on)
{
/*
	integrate equations of motion.
	
	INPUT:
	p		- class containing parameters.  following variablse
			  from that class used:

	display         - array telling which neurons are displayed.
	dt              - time step in ms.
	neurons         - number of neurons.
	no_display      - number of neurons that are displayed
	suffix		- suffix appended to output files.
	tmax            - maximum time for integrating, in ms.
	tmin            - starting time for integrating, in ms.

	INTEGRATION VARIABLES:
	z[0]    - phase
	z[1]	- excitatory drive (I_E)
	z[2]	- inhibitory drive (I_I)
*/
	// ---create derivs object
	derivs dzdt;

	// ---initialize dynamical variables.  this must be done now because
	//    it creates dzdt.no_vars, which is used elsewhere.
	float2 z;
	dzdt.init_z(p, z);

	// ---create runge_kutta object
	runge_kutta integrator(dzdt.no_vars, p.neurons);

	// ---reserve space
	int* who_spiked = (int*) calloc(p.neurons, sizeof(int));
	int* cnterr     = (int*) calloc(2, sizeof(int));

	// ---initialize spike count
	int* count = (int*) calloc(p.neurons, sizeof(int));
	int* count_total = (int*) calloc(3, sizeof(int));

	// ---initialize parametes in class dxdt;
	dzdt.set_params(p);

	// ---initialize time
	float t = p.tmin;

	int istep = 0;
	int* nspikes = (int*) calloc(p.no_mem+2, sizeof(int));
	int* aspikes = (int*) calloc(p.no_mem+2, sizeof(int));

	for (;;)
	{
		// ---exit if at end. include -dt because integrator
		//    increments t.
		if (t > p.tmax-p.dt) break;

		// --step equations of motion
		integrator.rk4(dzdt, z, t, p.dt);

		// ---modify parameters for any neuron that has spiked
		int n_spiked;
		spike(p, t, z, n_spiked, nspikes, who_spiked, cnterr);
		count_total[2] = p.no_mem ? nspikes[0] : 0;

		// ---write to output files and accumulate spike count
		if (t >= 0)
		{
			// --accumulate spikes
			for (int i=0; i < p.no_mem+2; i++)
				aspikes[i] += nspikes[i];

			// ---count spikes
			count_total[0]=count_total[1]=0;
			for (int i=0; i < n_spiked; i++)
			{
				count[who_spiked[i]]++;
				count_total[p.type[who_spiked[i]]]++;
			}

			// ---write
			if (write_on) write(p.display, n_spiked, p.no_display,
				t, who_spiked, count_total, aspikes, z, p);
		}
	}

	// ---compute average firing rate for each neuron
	float* nu = new float[p.neurons];
	for (int n=0; n < p.neurons; n++) nu[n] = count[n]/(p.tmax/1000);

	// ---print warnings about illegally crossing thresholds
	if (cnterr[0] > 0) write_warning("threshold at +pi illegally exceeded ",
		cnterr[0], " times; reduce time step!");
	if (cnterr[1] > 0) write_warning("threshold at -pi illegally exceeded ",
		cnterr[1], " times; reduce time step!");


	// ---free space
	free(who_spiked);
	free(cnterr);
	free(count);
	free(count_total);
	free(nspikes);

	return nu;
}

void spike(network& p, float t, float2& z, int& n_spiked, int* nspikes,
	int* who_spiked, int* cnterr)
{
/*
	perform various tasks if a neuron spiked

	INPUT:
	p		- class containing parameters
	t		- time in ms.
	z		- dynamical variables

	OUTPUT:
	n_spiked	- number of neurons that spiked on this time step
	who_spiked	- which neurons spiked on this time step
	cnterr		- cnterr[0] = number of times neurons went above +pi
			  cnterr[1] = number of times neurons went below -pi
	z		- gets updated if there is a spike
*/
	// ---constants
	static double pi = 4.0*atan(1.0);
	static double twopi = 8.0*atan(1.0);

	static int extra=1;
	static int n_dp=0, n_dp_nu, nnu_in, ext_drive=0;
	static double** dp = (double**) calloc(2, sizeof(double*));
	static double** dp_nu;
	static int* lnu_in = (int*) calloc(p.neurons, sizeof(int));
	if (n_dp == 0)
	{
		// ---spike count distribution for neurons being zapped
		n_dp=10000;
		for (int i=0; i < 2; i++) dp[i] = dpoisson(p.p_zap[i], n_dp);

		// ---check to see if there is external drive.
		for (int i=0; i < p.neurons; i++)
			ext_drive += p.nu_in[i] > 0 ? 1 : 0;

		if (ext_drive)
		{
			// ---spike count distribution for input firing rate
			n_dp_nu=1000;
			dp_nu = input_dist(p, lnu_in, nnu_in, n_dp_nu);
		}
	}

	// ---initialization
	n_spiked=0;
	for (int i=0; i < p.no_mem+2; i++) nspikes[i]=0;

	for (int j=0; j < p.neurons; j++) if (z.x[0][j] > pi)
	{
		// ---keep track of who spiked
		who_spiked[n_spiked++] = j;

		// ---reset phase
		z.x[0][j] -= twopi;

		// --update post-synaptic potentials
		update_post(p, z, nspikes, cnterr, j);
	}

	// ---extra neuron
	if (t > p.t_extra && extra)
	{
		// --update post-synaptic potentials
		update_post(p, z, nspikes, cnterr, p.j_extra);
		extra=0;
	}

	// ---bring any phase below -pi back into range.
	for (int j=0; j < p.neurons; j++)
	{
		if (z.x[0][j] < -pi)
		{
			cnterr[1]++;
			z.x[0][j] += twopi*((int) (-(z.x[0][j]-pi)/twopi));
		}
	}

	// ---external excitatory drive, if any
	if (ext_drive) for (int i=0; i < p.neurons; i++) z.x[1][i] +=
		p.wbar[0][p.type[i]]*dp_nu[lnu_in[i]][(int) (n_dp_nu*drand())];

	// ---zap: add pi to phase.
	if (p.no_mem && t >= p.t_zap_on && t < p.t_zap_off)
	{
		for (int i=0; i < p.num[0]; i++) if (p.who_zap[i])
			z.x[1][i] += p.wbar[0][0]*dp[0][(int) (n_dp*drand())];
	}

	// ---dezapping: multiply phase by small constant (p.dezap) to
	//    bring it near 0.
	if (p.no_mem && t >= p.t_dzap_on && t < p.t_dzap_off)
	{
		for (int i=0; i < p.num[0]; i++) if (p.who_zap[i])
			z.x[2][i] += p.wbar[1][0]*dp[1][(int) (n_dp*drand())];
	}

	if (t >= p.t_dzap_off && p.rezap)
	{
		p.rezap=0;
		float delta_zap = p.t_dzap_off - p.t_zap_on + p.dtzap;
		p.t_zap_on   += delta_zap;
		p.t_zap_off  += delta_zap;
		p.t_dzap_on  += delta_zap;
		p.t_dzap_off += delta_zap;
		for (int j=0; j < p.no_mem; j++)
			for (int k=0; k < p.num[0]; k++) p.who_zap[k] =
				p.eta[p.newzap][k] > 0 ? 1 : 0;
	}

}

void update_post(network& p, float2& z, int* nspikes, int* cnterr, int j)
{
	static double pi = 4.0*atan(1.0);
	static double twopi = 8.0*atan(1.0);

	// ---memories
	if (p.no_mem)
	{
		if (p.type[j] == 0) for (int i=0; i < p.nwmem[j]; i++)
			nspikes[p.wmem[j][i]]++;
	}

	// ---raw counts
	nspikes[p.no_mem + p.type[j]]++;

	// ---check to see if we have gone around an extra 2*pi in phase.
	if (z.x[0][j] > pi)
	{
		cnterr[0]++;
		z.x[0][j] -= twopi*((int) ((z.x[0][j]+pi)/twopi));
	}

	// ---loop over post-synaptic neurons
	int tpre=p.type[j];
	for (int m=0; m < p.wntot[j]; m++)
	{
		// ---update I_E or I_I, depending on type of
		//    presynaptic neuron.
		z.x[tpre+1][p.wn[j][m]] += p.w[j][m];
	}
}

double** input_dist(network& p, int* l, int& n, int n_dp)
{
/*
	INPUT:
	p.dt		- time step, in ms.
	p.dnu		- rates are quantized at dnu (nu_in = integer*dnu).
	p.nu_in[i]	- input firing rate of neuron i.
	n_dp_nu		- number of elements to use in approximation of
			  Poisson.

	OUTPUT:
	l[i]		- if l[i]=k, then firing rate = nu_k = integer*dnu,
			  i=0, ..., p.neurons-1.
	n		- number of elements in dp.

	RETURNS:
	dp_nu[l[i]][j]	- probability that dp_nu[l[i]][j]=n is
			  (nu_k dt)^n exp(-nu_k dt)/n!. j runs from 0 to n_dp-1.

*/
	// ---min and max
	float numin = min(p.neurons, p.nu_in);
	float numax = max(p.neurons, p.nu_in);

	// ---check
	if (numin < 0) write_err("numin < 0 in function input_dist.");

	// ---count number of unique intervals. first make a dummy
	//    variable big enough for all intervals.
	int ndum=(int) ((numax-numin)/p.dnu+1);
	int* ldum = (int*) calloc(ndum, sizeof(int));

	for (int i=0; i < p.neurons; i++)
		ldum[(int) ((p.nu_in[i] - numin)/p.dnu + 0.5)]=1;

	// ---count number of entries filled
	n=0;
	for (int i=0; i < ndum; i++) n += ldum[i];

	// ---make dp
	int cnt=0;
	double** dp = (double**) calloc(n, sizeof(double*));
	for (int i=0; i < ndum; i++) if (ldum[i])
		dp[cnt++]=dpoisson(p.dt*(p.dnu*i+numin)/1000, n_dp);

	// ---make l, which maps firing rate to index of dp. first
	//    step is to make ldum perform the mapping.
	cnt=0;
	for (int i=0; i < ndum; i++) ldum[i] = ldum[i] == 0 ? 0 : cnt++;

	// ---now, make l
	for (int i=0; i < p.neurons; i++)
		l[i]=ldum[(int) ((p.nu_in[i] - numin)/p.dnu + 0.5)];

	free(ldum);

	return dp;
}


void write(
	int* display, int n_spiked, int no_display, float t,
	int* who_spiked, int* count_total, int* aspikes, float2& z, network& p)
{
/*
	write output
*/
	// ---first step: open output files---

	// ---reserve space
	static int first = 1, n_outfile=7;
	static float twrite = 0;
	static FILE** f = (FILE**) calloc(n_outfile, sizeof(FILE*));
	static char** filename = cchar(n_outfile, 1024);

	// ---names of files
	sprintf(filename[0], "%s%s", "phase", p.suffix);
	sprintf(filename[1], "%s%s", "i_e", p.suffix);
	sprintf(filename[2], "%s%s", "i_i", p.suffix);
	sprintf(filename[3], "%s%s", "v", p.suffix);
	sprintf(filename[4], "%s%s", "s", p.suffix);
	sprintf(filename[5], "%s%s", "rate", p.suffix);
	sprintf(filename[6], "%s%s", "arate", p.suffix);


	// ---open files
	if (first)
	{
		first = 0;
		for (int i=0; i < n_outfile; i++)
			f[i] = fopen(filename[i], "w");
	}

	// ---second: output data---
	
	// ---i is important: it keeps track of which file we are on!
	int i;

	// ---phase and excitatory and inhibitory drive
	for (i=0; i < 3; i++)
	{
		fprintf(f[i], "%14.8f ", t/1000);
		for (int n=0; n < no_display; n++)
			fprintf(f[i], "%14.8f ", z.x[i][display[n]]);
		fprintf(f[i], "\n");
	}

	// ---voltage. no i++ because it was incremented at end of loop.
	//    choose epsilon so that vmax = 20.
	float vmax=20, vmin=-80;
	fprintf(f[i], "%14.8f ", t/1000);
	for (int n=0; n < no_display; n++)
	{
		int j = display[n];
		double v = p.vbar + p.deltav*tan(0.5*z.x[0][j]);
		v = v > vmax ? vmax : v;
		v = v < vmin ? vmin : v;
		fprintf(f[i], "%14.8f ", v);
	}
	fprintf(f[i], "\n");

	// ---spike times
	i++;
	fprintf(f[i], "%14.8f ", t/1000);
	for (int n=0; n < n_spiked; n++) fprintf(f[i], "%d ", who_spiked[n]);
	fprintf(f[i], "\n");

	// ---instantaneous rate. last quantity is (nu-f*nu_mem)/(1-f).
	i++;
	fprintf(f[i], "%14.8f %14.8f %14.8f %14.8f %14.8f\n", t/1000,
		count_total[0]/(p.num[0]*p.dt/1000),
		count_total[1]/(p.num[1]*p.dt/1000),
		count_total[2]/(p.f*p.num[0]*p.dt/1000),
		(count_total[0]-count_total[2])/((1-p.f)*p.num[0]*p.dt/1000));

	// ---average rate
	i++;
	double tnorm_e = p.num[0]*p.t_average/1000;
	double tnorm_i = p.num[1]*p.t_average/1000;
	if (t > twrite)
	{
		twrite += p.t_average;

		// ---excitatory and inhibitory rates
		fprintf(f[i], "%14.8f %14.8f %14.8f", t/1000,
			aspikes[p.no_mem]/tnorm_e, aspikes[p.no_mem+1]/tnorm_i);

		// ---memories
		for (int j=0; j < p.no_mem; j++)
			fprintf(f[i], " %14.8f", aspikes[j]/(p.f*tnorm_e));

		// ---CR
		fprintf(f[i], "\n");

		// ---zero out spikes
		for (int i=0; i < p.no_mem+2; i++) aspikes[i]=0;
	}
}

void count_ea(network& p) 
{
/*
	count number of endogenously active cells and print results
*/
	int n_ea[2] = {0, 0};
	for (int i=0; i < p.neurons; i++) if (p.mu[i] > 0.0) n_ea[p.type[i]]++;

	fprintf(stderr, "%5i out of %6i %s\n", n_ea[0], p.num[0],
		"excitatory neurons endogenously active");
	fprintf(stderr, "%5i out of %6i %s\n", n_ea[1], p.num[1],
		"inhibitory neurons endogenously active");
}

void write_filelist()
{
	write_message("theta1f.TMP_MORE_TMP",
"\n"
"output files (followed by '.suffix', if suffix is on command line).\n"
"'e and i' means file_e contains excitatory cells, file_i contains inhibitory.\n"
"\n"
"arate	- population averaged rate (Hz), averaged over t_average ms, vs time.\n"
"	  columns are: time (s), nu_E, nu_I, nu_j, j=1, no_mem.\n"
"eta	- value of eta for each memory. not printed if no memories.\n"
"f	- population averaged rates (Hz) for excitatory and inhibitory cells.\n"
"fc	- firing rate (Hz) of each neuron. e and i.\n"
"i	- drive versus time for displayed neurons. e and i.\n"
"mu	- external drive for each cell: current (mV) and mu. e and i.\n"
"phase	- phase versus time for displayed neurons.\n"
"rate	- population averaged rate (Hz), averaged over 1 bin, vs time.\n"
"	  columns are: time (s), nu_E, nu_I, nu_mem, (nu_E - f*nu_mem)/(1-f).\n"
"s	- list of who spiked on each time step.\n"
"v	- voltage (mV) versus time for displayed neurons.\n"
"\n");

exit(1);
}

void write_comline()
{
	write_message("theta1d.TMP_MORE_TMP",
"\n"
"command line:  theta1d input [suffix] -flags\n"
"\n"
"  Network of theta-neurons.\n"
"  \n"
"  If suffix is present (must be third argument, and not one of the\n"
"  flags), '.suffix' will be appended to all output files.\n"
"\n"
"  ----flags\n"
"  -h: this message.\n"
"  -i: list of input parameters.\n"
"  -c: description of how current distribution is computed.\n"
"  -s: get a sample input file that works pretty well.\n"
"  -o: list of (generally dimensionless) parameters generated by the code.\n"
"  -f: list of output files.\n"
"  -w: if nonzero, write weights to files w_n where n is presynaptic neuron\n"
"      (as specified in display). format is: n_post weight. also make\n"
"      histograms with arg(-w) bins and write to files w_e and w_i. format is:\n"
"      weight count [count]. w_i has one count. w_e has two: weights from\n"
"      non-memory and memory neurons, respectively.\n"
"\n");

exit(1);
}
