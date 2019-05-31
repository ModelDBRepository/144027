#include "params.h"
#include <cstring>
/*
	Peter Latham
        September, 2005
	version 1
*/

// constructor
params::params(int brgc, char** brgv, int prnt_cur, int prnt_in, int prnt_s)
{
	// ---message explaining variables.
	if (prnt_in) write_comline1();

	// ---sample list of variables.
	if (prnt_s) write_comline2();

	// ---message explaining how current distribution is constructed.
	if (prnt_cur) write_comline3();

	// ---get elements of input file.
	charplus z = get_chars(fopen(brgv[1], "r"));

	// ---check
	if (z.n == 0) write_err("file \"", brgv[1], "\" does not exist.");

	// ---count elements
	int argc=0;
	for (int i=0; i < z.n; i++) argc += z.m[i];

	// ---reserve space
	char** argv = (char**) calloc(argc, sizeof(char*));

	// ---transfer data
	argc=0;
	for (int i=0; i < z.n; i++) for (int j=0; j < z.m[i]; j++)
	{
		argv[argc] = (char*) calloc(strlen(z.x[i][j])+1, sizeof(char));
		strcpy(argv[argc], z.x[i][j]);
		argc++;
	}

	// ------read variables------
	int nmatch;
	float dum=1;
	float pi = 4.0*atan(1.0);

	// ---miscellaneous---
	//    note: no_display needs to be positive for call to parse_comline.
	frac_inh    = parse_comline(argc, argv, "-frac_inh", 1, (float) 0.2);
	neurons     = parse_comline(argc, argv, "-neurons", 1, 1000);
	seed        = parse_comline(argc, argv, "-seed",    1, 1);
	axons[0][0] = parse_comline(argc, argv, "-axons_e", 1, (float) 1000);
	axons[0][1] = parse_comline(argc, argv, "-axons_e", 2, (float) 1000);
	axons[1][0] = parse_comline(argc, argv, "-axons_i", 1, (float) 1000);
	axons[1][1] = parse_comline(argc, argv, "-axons_i", 2, (float) 1000);
	display     = parse_comline(argc, argv, "-display", no_display=1);

	// ---basic time parameters, including time constants---
	dt        = parse_comline(argc, argv, "-dt",        1, (float) 0);
	tmin      = parse_comline(argc, argv, "-tmin",      1, (float) 0)*1000;
	tmax      = parse_comline(argc, argv, "-tmax",      1, (float) 1)*1000;
	tau_m[0]  = parse_comline(argc, argv, "-tau",       1, (float) 10);
	tau_m[1]  = parse_comline(argc, argv, "-tau",       2, (float) 10);
	tau_s[0]  = parse_comline(argc, argv, "-tau_s",     1, (float) 5);
	tau_s[1]  = parse_comline(argc, argv, "-tau_s",     2, (float) 5);
	t_average = parse_comline(argc, argv, "-t_average", 1, (float) 100);

	// ---voltages---
	e_exc = parse_comline(argc, argv, "-ureverse", 1, (float)  0.0);
	e_inh = parse_comline(argc, argv, "-ureverse", 2, (float) -80.0);
	v_r   = parse_comline(argc, argv, "-v_r",      1, (float) -65.0);
	v_t   = parse_comline(argc, argv, "-v_t",      1, (float) -50.0);

	// ---current distribution---
	ia_mean = parse_comline(argc, argv, "-ia_mean", dum, ia_no, "");
	ia_sd   = parse_comline(argc, argv, "-ia_sd",   dum, nmatch, "");
	ia_min  = parse_comline(argc, argv, "-ia_min",  dum, nmatch, "");
	ia_max  = parse_comline(argc, argv, "-ia_max",  dum, nmatch, "");
	ia_p    = parse_comline(argc, argv, "-ia_p",    dum, nmatch, "");

	// ---input firing rate distribution---
	nu_mean = parse_comline(argc, argv, "-nu_mean", dum, nu_no, "");
	nu_sd   = parse_comline(argc, argv, "-nu_sd",   dum, nmatch, "");
	nu_min  = parse_comline(argc, argv, "-nu_min",  dum, nmatch, "");
	nu_max  = parse_comline(argc, argv, "-nu_max",  dum, nmatch, "");
	nu_p    = parse_comline(argc, argv, "-nu_p",    dum, nmatch, "");
	dnu     = parse_comline(argc, argv, "-dnu",     1, (float) 10);

	// ---random component of weights---
	dw          = parse_comline(argc, argv, "-dw",     1, (float) 0);
	vpsp[0][0] = parse_comline(argc, argv, "-epsp", 1, (float) 1.0);
	vpsp[0][1] = parse_comline(argc, argv, "-epsp", 2, (float) 1.0);
	vpsp[1][0] = parse_comline(argc, argv, "-ipsp", 1, (float) -1.5);
	vpsp[1][1] = parse_comline(argc, argv, "-ipsp", 2, (float) -1.5);

	// ---memories---
	zap        = parse_comline(argc, argv, "-zap",      1, (float) pi);
	dzap       = parse_comline(argc, argv, "-dzap",     1, (float) 0.1);
	eps        = parse_comline(argc, argv, "-eps",      1, (float) 0.1);
	pspmin     = parse_comline(argc, argv, "-psprange", 1, (float) 0.0);
	pspmax     = parse_comline(argc, argv, "-psprange", 2, (float) 3.0);
	f          = parse_comline(argc, argv, "-f",        1, (float) 0.1);
	norm       = parse_comline(argc, argv, "-norm",     1, 0);
	whozap     = parse_comline(argc, argv, "-whozap",   1, 1);
	no_mem     = parse_comline(argc, argv, "-no_mem",   1, 1);
	nu_zap[0]  = parse_comline(argc, argv, "-nu_zap",   1, (float) 10);
	nu_zap[1]  = parse_comline(argc, argv, "-nu_zap",   2, (float) -1);
	nu_zap[1]  = nu_zap[1] == -1 ? nu_zap[0] : nu_zap[1];
	t_zap_on   = parse_comline(argc, argv, "-t_zap",    1, (float) 0)*1000;
	t_zap_off  = parse_comline(argc, argv, "-t_zap",    2, (float) 0)*1000;
	t_dzap_on  = parse_comline(argc, argv, "-t_dzap",   1, (float) 0)*1000;
	t_dzap_off = parse_comline(argc, argv, "-t_dzap",   2, (float) 0)*1000;
	dneps      = parse_comline(argc, argv, "-deps", (float) 1, no_deps, "");

	// ---zap a second time. kind of a hack.
	rezap      = parse_comline(argc, argv, "-rezap",    0, 0);
	dtzap      = parse_comline(argc, argv, "-dtzap",    1, (float) 100);
	newzap     = parse_comline(argc, argv, "-newzap",   1, 2);

	// ---extra spike---
	t_extra    = parse_comline(argc, argv, "-t_extra", 1, (float) 6)*1000;
	j_extra    = parse_comline(argc, argv, "-j_extra", 1, 0);


	// ---output files---

	// ---get suffix
	get_suffix(brgc, brgv);

	// ---some checking---

	// ---no of neurons displayed
	if (no_display > neurons)
	{
		write_warning("number of displayed neurons reduced to ",
			neurons, ".");
		no_display = neurons;
	}

	// ---which neurons displayed
	for (int i=0; i < no_display; i++) if (display[i] >= neurons)
	{
		write_warning("display neuron #", i+1, " set to zero.");
		display[i]=0;
	}

	// ---who gets zapped
	if (no_mem && whozap > no_mem)
	{
		write_warning("whozap reduced to ", no_mem, ".");
		whozap=no_mem;
	}

	// ---psprange
	if (pspmin < 0)
	{
		write_warning("pspmin set to 0");
		pspmin=0;
	}

	// ---reduce whoazap and newzap so that 0 refers to memory 1, etc.
	//    negative means zap random neurons.
	whozap--;
	newzap--;
}

void params::get_suffix(int brgc, char** brgv)
{
/*
	get suffix.
	if brgc > 2, set suffix to ".brgv[2]". otherwise, set it to "".
	also, if brgv[2][0] = '-', assume it is a flag and set suffix to "".
*/
	suffix = (char*) calloc(1024, sizeof(char));
	if (brgc > 2)
	{
		if (brgv[2][0] != '-') sprintf(suffix, ".%s", brgv[2]);
		else sprintf(suffix, "");
	}
}

// ---message explaining variables.
void params::write_comline1()
{
	write_message("theta11.TMP_MORE_TMP",
"\n"
"input parameters. parameters followed by [i] require two arguments.\n"
"the first is excitatory, the second inhibitory.\n"
"\n"
"---miscellaneous\n"
"frac_inh	- fraction of neurons that are inhibitory.\n"
"neurons		- number of neurons.\n"
"axons_e[i]	- mean number of axons from exc. neurons to neurons of type i.\n"
"axons_i[i]	- mean number of axons from inh. neurons to neurons of type i.\n"
"display[i]	- array telling which neurons are displayed.\n"
"seed		- seed for random number generator.\n"
"\n"
"---basic time parameters, including time constants.\n"
"dt		- time step in ms.\n"
"tmin		- starting time of simulation, in s.\n"
"tmax		- maximum time of simulation, in s.\n"
"tau[i]		- membrane time constant, ms.\n"
"tau_s[i]	- time constant for excitatory synapses, ms.\n"
"t_average	- time to average over for output, ms.\n"
"\n"
"---voltages\n"
"ureverse[i]	- excitatory reversal potential, mV.\n"
"v_r		- nominal resting membrane potential, mV.\n"
"v_t		- nomianl threshold, mV.\n"
"\n"
"---current distribution (type \"theta1 -c\" for details)\n"
"\n"
"ia_mean[i]	- mean of each Gaussian distribution.\n"
"ia_sd[i]	- standard deviation of each distribution before\n"
"		  truncation (see following two parameters).\n"
"ia_min[i]	- minimum value of each distribution.\n"
"ia_max[i]	- maximum value of each distribution.\n"
"ia_p[i]		- probability that a neuron will come from a particular\n"
"		  Gaussian. positive for excitatory distribution, negative\n"
"		  for inhibitory.\n"
"\n"
"---input firing rate distribution (type \"theta1 -c\" for details)\n"
"\n"
"nu_mean[i]	- mean of each Gaussian distribution.\n"
"nu_sd[i]	- standard deviation of each distribution before\n"
"		  truncation (see following two parameters).\n"
"nu_min[i]	- minimum value of each distribution.\n"
"nu_max[i]	- maximum value of each distribution.\n"
"nu_p[i]		- probability that a neuron will come from a particular\n"
"		  Gaussian. positive for excitatory distribution, negative\n"
"		  for inhibitory.\n"
"\n"
"---random component of weights\n"
"dw		- each connection strength is multiplied by (1+xi)\n"
"		  where xi is a random variable uniformly distributed\n"
"		  between -dw and dw.\n"
"epsp[i]		- EPSP for post-synaptic cell of type i.\n"
"ipsp[i]		- IPSP for post-synaptic cell of type i.\n"
"\n"
"\n"
"---memories\n"
"no_mem		- number of memories.\n"
"zap		- add phase zap to any neuron that is zapped.\n"
"dzap		- phase gets mulitplied by dezap, which must be less than 1.\n"
"		  this brings phase near zero and turns neuron off.\n"
"eps		- strength of memory.\n"
"psprange	- range of allowed EPSPs, in mV. two args: min and max.\n"
"f		- fraction of neurons involved in a memory.\n"
"norm		- 0: both pre and post-synaptic normalization,\n"
"		  1: postsynatpic normalization only.\n"
"whozap		- labels which memory getz zapped. 1 is first memory,\n"
"		  2 is second, etc. 0 or negative means zap random neurons.\n"
"nu_zap[i]	- firing rate of on (i=0) and off (i=1) zapping neurons, Hz.\n"
"t_zap_on	- time in s when zapping starts.\n"
"t_zap_off	- time in s when zapping stops.\n"
"t_dzap_on	- time in s when dezapping starts.\n"
"t_dzap_off	- time in s when dezapping stops.\n"
"deps		- memories are multiplied by deps. format is (m deps) where\n"
"		  m is the memory and deps is the multiplier. this can be\n"
"		  repeated as many times as you want, and it doesn't have to\n"
"		  be all on one line. first memory is 1, not 0.\n"
"\n"
"---zap a second time\n"
"rezap		- if on, zap the neurons a second time and use the following:\n"
"dtzap		- zapping sequence starts over after time deltay dtzap, in ms.\n"
"newzap		- new memory to zap.\n"
"\n"
"---extra spike\n"
"t_extra		- time at which extra spike appears, in s.\n"
"j_extra		- neuron on which extra spike appears.\n"
"\n");

exit(1);
}

// ---sample list of variables.
void params::write_comline2()
{
	write_message("theta12.TMP_MORE_TMP",
"\n"
"# ---miscellaneous\n"
"-frac_inh	0.20\n"
"-neurons	1000\n"
"-axons_e	100 100\n"
"-axons_i	100 100\n"
"-display	1 7 19 23\n"
"-seed		4\n"
"\n"
"# ---basic time parameters, including time constants.\n"
"-dt		1\n"
"-tmin		0\n"
"-tmax		10\n"
"-tau		10 10\n"
"-tau_s		5 5\n"
"-t_average	100\n"
"\n"
"# ---voltages\n"
"-ureverse	0 -70\n"
"-v_r		-65\n"
"-v_t		-50\n"
"\n"
"# ---current distribution\n"
"-ia_mean	3.0\n"
"-ia_sd		100\n"
"-ia_min		0\n"
"-ia_max		4.7\n"
"-ia_p		1\n"
"\n"
"# ---input firing rate distribution\n"
"-nu_mean	100\n"
"-nu_sd		50\n"
"-nu_min		0\n"
"-nu_max		200\n"
"-nu_p		1\n"
"\n"
"# ---random component of weights\n"
"-dw		0.00\n"
"-epsp		 1.0  1.0\n"
"-ipsp		-1.5 -1.5\n"
"\n"
"# ---memories\n"
"-no_mem		10\n"
"-zap		3.14159\n"
"-dzap		0.1\n"
"-eps		0.1\n"
"-psprange	0.0 3.0\n"
"-f		0.1\n"
"-norm		1\n"
"-whozap		0\n"
"-nu_zap		10 10\n"
"-t_zap		2 2.1\n"
"-t_dzap		3 3.1\n"
"-deps\n"
"1 1.1\n"
"3 0.9\n"
"\n"
"# ---extra spike\n"
"-t_extra	6.0\n"
"-j_extra	0\n"
"\n");

exit(1);
}

// ---message explaining how current distribution is constructed.
void params::write_comline3()
{
	write_message("theta13.TMP_MORE_TMP",
"\n"
"The distribution of currents or firing rates is the weighted sum of\n"
"truncated Gaussian distributions. The variable that determines the\n"
"weight of each distribution is p (either ia_p, for current, or nu_p,\n"
"for firing rate), as follows:\n"
"\n"
"If all the p are non-negative, then the excitatory and inhibitory\n"
"neurons have the same distribution.\n"
"\n"
"If any of the p are negative, then the excitatory neurons go with\n"
"the distributions with non-negative p and the inhibitory neurons go\n"
"with the distributions with negative p. For both types the total\n"
"probability is forced to add to 1, and of course |p| will be used\n"
"to compute the probability.  \n"
"\n");

exit(1);
}
