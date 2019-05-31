#include "lib.h"
#include <cstring>
/*
	new --> calloc and delete --> free:
	May 20, 2004: 
		int** get_cols

	free --> delete:
	May 20, 2004:
		zfree

	ttd:
	1.	in fgets, check for overrun.

--------------------------------------------------------------------------------

	cat2		- cat together two strings.
	choose(k,n)	- n choose k.
	drand		- random number generator.
	entropy		- naive estimate of entropy.
	entropy2	- naive estimate of second moment of log(p).
	get_chars	- parse file and return list of character
			  strings pointing to each element on each line.
	get_cols	- get columns from file (including stdin).
	header_info	- get info from headers in ".dat" files.
	isnumber	- turns string into a number with lots of checking.
	kl_divergence	- Kullback-Leibler distance between two distributions.
	length_string	- returns length of \0 terminated string, including \0.
	max		- finds maximum of array.
	min		- finds minimum of array.
	mmult		- multiply matrix times matrix or matrix times vector.
	cint		- 2-D integer array, using calloc.
	cfloat		- 2-D float array, using calloc.
	cdouble		- 2-D double array, using calloc.
	newchar		- 2-D character array.
	newdouble	- 2-D double array, using new.
	newfloat	- 2-D float array.
	newint		- 2-D integer array.
	init_constant	- initizlize 1 and 2 D arrays.
	mi		- computes mutual information
	rabs		- absolute value.
	read_binary	- read binary file
	read_n		- read n lines of a file and return data in kth column.
	read_plexon	- read plexon file.
	readline	- read line of file.
	parse_c		- parse data file -- return pointer to string.
	parse_comline	- extract info from command line.
	s_parse		- simple string parse.
	unit_no		- return cell number in unit format given
			  input in _either_ plexon or unit format.
	write_err	- write an error message and exit(1).
*/

char* cat2(char* prefix, char* suffix)
{
/*
	cat together two strings and return the result
*/
	char* dum = new char[80];
	strcpy(dum, prefix);
	strcat(dum, suffix);
	return dum;
}

int choose(int k, int n)
{
/*
	compute n!/(k! (n-k)!)

	will get overflow if n and k are too big.
*/
	if (k > n) return 0;

 	int c = 1;
	k = k > n/2 ? n-k: k;
	for (int i=0; i < k; i++) c *= (n-i);
	for (int i=0; i < k; i++) c /= (i+1);

	return c;
}

double drand()
{
	// return a random number between 0 and 1
	// use srand(seed) to start seed.
	// 4/4/00:  RAND_MAX --> RAND_MAX-1, so 1 cannot be returned.
	// 8/4/04:  RAND_MAX --> RAND_MAX. why I thought RAND_MAX-1
	//                       would fix this is a mystery. would
	//                       like to go to RAND_MAX+1, but this
	//                       causes overflow and returns negative
	//                       numbers.
	return (double) rand()/RAND_MAX;
}

double drand(float x)
{
	// return a random number between 0 and x
	return x*rand()/RAND_MAX;
}

double drand(float xmin, float xmax)
{
	// return a random number between xmin and xmax
	return xmin + (xmax-xmin)*rand()/RAND_MAX;
}

double entropy(int n, int* p)
{
/*
	overloaded function (int* p --> double s)

	compute entropy of a probability distribution, p.
*/
	// ---find ntot
	float ntot=0;
	for (int i=0; i < n; i++) ntot += p[i] > 0 ? p[i] : 0;

	double s=0;
	for (int i=0; i < n; i++)
		s -= p[i] > 0 ? (p[i]/ntot)*log(p[i]/ntot) : 0;

	return s/log(2);
}

double entropy(int n, float* p, int norm)
{
/*
	overloaded function (float* p --> double s)

	compute entropy of a probability distribution, p.
*/
	double s=0;
	for (int i=0; i < n; i++) s -= p[i] > 0 ? p[i]*log(p[i]) : 0;

	if (norm)
	{
		// normalize probability distribution
		double z=0;
		for (int i=0; i < n; i++) z += p[i] > 0 ? p[i] : 0;
		if (z > 0) s = s/z + log(z);
	}

	return s/log(2);
}

double entropy(int n, double* p, int norm)
{
/*
	overloaded function (double* p --> double s)

	compute entropy of a probability distribution, p.
*/
	double s=0;
	for (int i=0; i < n; i++) s -= p[i] > 0 ? p[i]*log(p[i]) : 0;

	if (norm)
	{
		// normalize probability distribution
		double z=0;
		for (int i=0; i < n; i++) z += p[i] > 0 ? p[i] : 0;
		if (z > 0) s = s/z + log(z);
	}

	return s/log(2);
}

double entropy2(int n, int* p)
{
/*
	overloaded function (int* p --> double s)

	compute second moment of log p.
*/
	// ---find ntot
	float ntot=0;
	for (int i=0; i < n; i++) ntot += p[i] > 0 ? p[i] : 0;

	double s=0;
	for (int i=0; i < n; i++)
		s += p[i] > 0 ? (p[i]/ntot)*pow(log(p[i]/ntot), 2) : 0;

	return s/(log(2)*log(2));
}

double entropy2(int n, float* p, int norm=0)
{
/*
	overloaded function (float* p --> double s)

	compute second moment of log p.
*/
	double z=1;
	if (norm)
	{
		// normalize probability distribution
		double z=0;
		for (int i=0; i < n; i++) z += p[i] > 0 ? p[i] : 0;
	}

	double s=0;
	for (int i=0; i < n; i++)
		s += p[i] > 0 ? (p[i]/z)*pow(log(p[i]/z), 2) : 0;

	return s/(log(2)*log(2));
}

double entropy2(int n, double* p, int norm=0)
{
/*
	overloaded function (double* p --> double s)

	compute second moment of log p.
*/
	double z=1;
	if (norm)
	{
		// normalize probability distribution
		double z=0;
		for (int i=0; i < n; i++) z += p[i] > 0 ? p[i] : 0;
	}

	double s=0;
	for (int i=0; i < n; i++)
		s += p[i] > 0 ? (p[i]/z)*pow(log(p[i]/z), 2) : 0;

	return s/(log(2)*log(2));
}

charplus get_chars(FILE* f)
{
/*
	read from file f, except skips lines whose first element is #.
	delimiter is white space or tab.

	RETURNS:
	charplus z:
	z.x[i][j]	- jth element in line i, in char* format
	z.m[i]		- number of elements on line i
	z.n		- number of lines.  returns 0 if the file
			  does not exist.
*/
	// ---declarations
	struct charplus z;

	// ---check that file exists.  if not, return with z.n = 0.
	//    reserve space just in case.
	if (!f)
	{
		z.n=0;
		z.m = new int[1];
		z.x = new char**[1];
		z.x[0] = new char*[1];
		z.x[0][0] = new char[1];
		return z;
	}

	// reserve space and make a few declarations
	char* line = new char[MLL_lib];
	int n_elements;
	ostrstream m_ptr;

	strstream inout;
	int line_no=0;
	while (fgets(line, MLL_lib, f))
	{
		inout << line;
		line_no++;
	}

	// ---reserve slightly too much memory, since some lines are
	//    comments.
	z.x = new char**[line_no];
	z.n = 0;
	for (int i=0; i < line_no; i++)
	{
		inout.getline(line, MLL_lib);

		// parse
		int* nw = parse_c(line, n_elements);

		// check that there is at least one element on the line
		if (n_elements == 0)
		{
			z.x[z.n] = new char*[1];
			z.x[z.n][0] = new char[1];
			z.n++;
			m_ptr.write((char*) &n_elements, sizeof(n_elements));
		}

		// check for comment
		else if (line[nw[0]] != '#')
		{
			// ---data exists

			// reserve space
			z.x[z.n] = new char*[n_elements];
			for (int i=0; i < n_elements; i++)
			{
				z.x[z.n][i] =
					new char[length_string(&line[nw[i]])];
				strcpy(z.x[z.n][i], &line[nw[i]]);
			}

			m_ptr.write((char*) &n_elements, sizeof(n_elements));
			z.n++;
		}

		delete [] nw;
	}

	delete [] line;

	z.m = (int*) m_ptr.str();
	return z;
}

float** get_cols(int n_col, int* col, int& n, FILE* f, int silent)
{
/*
	read from file f.  skips lines whose first element is #

	INPUT:
	n_col	- see col below.
	col	- reads columns col[0], col[1], ..., col[n_col-1].
		  skips lines with fewer than max(col) elements.
	f	- file.  in calling program:
		  FILE *f;
		  f = fopen("filename", "r");
		  or
		  f = stdin;
	silent	- if nonzero, do not announce empty lines


	OUTPUT:
	n	- number of columns

	RETURNS:
	x	- x[n][j] = jth element in column col[n]
*/

	// some checking
	int err=0;
	if (n_col < 1)
	{
		cerr << "function get_cols:  first argument is "
			<< "less than 1.  code will probably bomb." << endl;
		err=1;
	}
	for (int i=0; i < n_col; i++)
	{
		if (col[i] < 0)
		{
			cerr << "function get_cols:  col[" << i << "]"
				<< " is less than 0.  code will probably bomb."
				<< endl;
			err=1;
		}
	}
	if (err)
	{
		float** x = new float*[1];
		x[0] = new float[1];
		return x;
	}

	// reserve space
	float** x = new float*[n_col];

	// find maximum column
	int colmax = col[0];
	for (int i=1; i < n_col; i++)
		colmax = col[i] > colmax ?  col[i] : colmax;
	colmax += 1;

	// reserve space and make a few declarations
	char* line = new char[MLL_lib];
	int line_no=0, n_elements;
	ostrstream* x_ptr = new ostrstream[n_col];

	n=0;
	while (fgets(line, MLL_lib, f))
	{
		line_no++;

		// parse
		int* nw = parse_c(line, n_elements);

		// check for empty line
		if (n_elements == 0) {}

		// check for comment
		else if (line[nw[0]] == '#') {}

		// check for enough elements.
		else if (n_elements < colmax)
		{
			if (!silent)
				cerr << "line " << line_no << ":  fewer than "
					<< colmax << " elements; skipping"
					<< endl;
		}
		else
		{
			n++;

			for (int i=0; i < n_col; i++)
			{
				float y = atof(&line[nw[col[i]]]);
				x_ptr[i].write((char*) &y, sizeof(y));
			}
		}

		delete [] nw;
	}

	for (int i=0; i < n_col; i++) x[i] = (float*) x_ptr[i].str();

	delete [] line;
	fclose(f);

	return x;
}

int** get_cols(
	FILE* f, int& lines, int& n_elements, int flag, int min_elements)
{
/*
	read from file f, except skips lines whose first element is #.

	INPUT:
	f		- file to be read from.
	flag		- determines orientation of array x; see below.

	OPTIONAL INPUT:
	min_elements	- only lines with at least min_elements will
			  be saved. default=0

	OUTPUT:
	lines		- number of lines
	n_elements	- number of elements on each line.

	RETURNS:
	x[i][j]		- flag=1:  element i, line j.
			  flag=0:  line i, element j.
*/
	int** x;

	// ---check that file exists.  if not, return with lines =
	//    n_elements=0;
	if (!f)
	{
		lines=n_elements=0;
		return x;
	}

	// ---dummy variable; not used.
	char** dat;

	ostrstream x_ptr;
	ostrstream m_ptr;
	parse_file(f, x_ptr, m_ptr, lines, 0, dat, -1.0, ' ', min_elements);

	// ---check
	if (lines == 0)
	{
		n_elements=0;
		return x;
	}

	// ---number of elements on each line
	int* m = (int*) m_ptr.str();

	// ---check that m[i] is the same for all lines.  if it isn't,
	//    return with lines=-1.
	n_elements = m[0];
	for (int i=1; i < lines; i++) if (m[i] != n_elements)
	{
		lines=-1;
		return x;
	}

	// ---temporarily store the whole file
	double* xtmp = (double*) x_ptr.str();

	// ---transfer to x

	int cnt=0;
	if (flag)
	{
		// ---x[i][j] is element i, line j
		x = cint(n_elements, lines);

		for (int i=0; i < lines; i++) for (int j=0; j < n_elements; j++)
			x[j][i] = (int) floor(xtmp[cnt++] + 0.5);
	}
	else
	{
		// ---x[i][j] is line i, element j
		x = cint(lines, n_elements);

		for (int i=0; i < lines; i++) for (int j=0; j < n_elements; j++)
			x[i][j] = (int) floor(xtmp[cnt++] + 0.5);
	}

	delete [] xtmp;
	fclose(f);

	return x;
}

float** get_cols(
	FILE* f, int& lines, int& n_elements, float flag, int min_elements)
{
/*
	read from file f, except skips lines whose first element is #.

	INPUT:
	f		- file to be read from.
	flag		- determines orientation of array x; see below.

	OPTIONAL INPUT:
	min_elements	- only lines with at least min_elements will
			  be saved. default=0

	OUTPUT:
	lines		- number of lines
	n_elements	- number of elements on each line.

	RETURNS:
	x[i][j]		- flag=1:  element i, line j.
			  flag=0:  line i, element j.
*/
	float** x;

	// ---check that file exists.  if not, return with lines =
	//    n_elements=0;
	if (!f)
	{
		lines=n_elements=0;
		return x;
	}

	// ---dummy variable; not used.
	char** dat;

	ostrstream x_ptr;
	ostrstream m_ptr;
	parse_file(f, x_ptr, m_ptr, lines, 0, dat, -1.0, ' ', min_elements);

	// ---check
	if (lines == 0)
	{
		n_elements=0;
		return x;
	}

	// ---number of elements on each line
	int* m = (int*) m_ptr.str();

	// ---check that m[i] is the same for all lines.  if it isn't,
	//    return with lines=-1.
	n_elements = m[0];
	for (int i=1; i < lines; i++) if (m[i] != n_elements)
	{
		lines=-1;
		return x;
	}

	// ---temporarily store the whole file
	double* xtmp = (double*) x_ptr.str();

	// ---transfer to x

	int cnt=0;
	if (flag)
	{
		// ---x[i][j] is element i, line j
		x = newfloat(n_elements, lines);

		for (int i=0; i < lines; i++) for (int j=0; j < n_elements; j++)
			x[j][i] = xtmp[cnt++];
	}
	else
	{
		// ---x[i][j] is line i, element j
		x = newfloat(lines, n_elements);

		for (int i=0; i < lines; i++) for (int j=0; j < n_elements; j++)
			x[i][j] = xtmp[cnt++];
	}

	delete [] xtmp;
	fclose(f);

	return x;
}

double** get_cols(
	FILE* f, int& lines, int& n_elements, double flag, int min_elements)
{
/*
	read from file f, except skips lines whose first element is #.

	INPUT:
	f		- file to be read from.
	flag		- determines orientation of array x; see below.

	OPTIONAL INPUT:
	min_elements	- only lines with at least min_elements will
			  be saved. default=0

	OUTPUT:
	lines		- number of lines
	n_elements	- number of elements on each line.

	RETURNS:
	x[i][j]		- flag=1:  element i, line j.
			  flag=0:  line i, element j.
*/
	double** x;

	// ---check that file exists.  if not, return with lines =
	//    n_elements=0;
	if (!f)
	{
		lines=n_elements=0;
		x = (double**) calloc(1, sizeof(double*));
		x[0] = (double*) calloc(1, sizeof(double));
		return x;
	}

	// ---dummy variable; not used.
	char** dat;

	ostrstream x_ptr;
	ostrstream m_ptr;
	parse_file(f, x_ptr, m_ptr, lines, 0, dat, -1.0, ' ', min_elements);

	// ---check
	if (lines == 0)
	{
		n_elements=0;
		return x;
	}

	// ---number of elements on each line
	int* m = (int*) m_ptr.str();

	// ---check that m[i] is the same for all lines.  if it isn't,
	//    return with lines=-1.
	n_elements = m[0];
	for (int i=1; i < lines; i++) if (m[i] != n_elements)
	{
		lines=-1;
		return x;
	}

	// ---temporarily store the whole file
	double* xtmp = (double*) x_ptr.str();

	// ---transfer to x

	int cnt=0;
	if (flag)
	{
		// ---x[i][j] is element i, line j
		x = newdouble(n_elements, lines);

		for (int i=0; i < lines; i++) for (int j=0; j < n_elements; j++)
			x[j][i] = xtmp[cnt++];
	}
	else
	{
		// ---x[i][j] is line i, element j
		x = newdouble(lines, n_elements);

		for (int i=0; i < lines; i++) for (int j=0; j < n_elements; j++)
			x[i][j] = xtmp[cnt++];
	}

	delete [] xtmp;
	fclose(f);

	return x;
}

floatplus get_cols(FILE* f)
{
/*
	read from file f, except skips lines whose first element is #.
	call:
		floatplus z = get_cols(stdin);
	or (I believe)
		floatplus z = get_cols(fopen("filename", "r"));

	INPUT:
	f	- input file.

	RETURNS:
	floatplus z:
	z.x[i][j]	- jth element in line i
	z.m[i]		- number of elements on line i
	z.n		- number of lines.  returns 0 if the file
			  does not exist.
*/
	// ---declarations
	struct floatplus z;

	// ---check that file exists.  if not, return with z.n = 0.
	//    reserve space just in case.
	if (!f)
	{
		z.n=0;
		z.m = new int[1];
		z.x = new double*[1];
		z.x[0] = new double[1];
		z.ndat=1;
		z.ddat = new double[1];
		return z;
	}

	// reserve space and make a few declarations
	ostrstream x_ptr;
	ostrstream m_ptr;
	char** dat;

	z.ddat = parse_file(f, x_ptr, m_ptr, z.n, 0, dat, 0.0, ' ', 0);

	z.m = (int*) m_ptr.str();

	// ---temporarily store the whole file
	double* x = (double*) x_ptr.str();

	// ---transfer to z.x
	z.x = newdouble(z.n, z.m);

	int cnt=0;
	for (int i=0; i < z.n; i++) for (int j=0; j < z.m[i]; j++)
		z.x[i][j] = x[cnt++];

	delete [] x;
	fclose(f);

	return z;
}

floatplus get_cols(
	FILE* f, int ndat, char** dat, double dflt=-1.0, char sep=' ',
	int min_elements)
{
/*
	read from file f, except skips lines whose first element is #.
	call:
		floatplus z = get_cols(stdin);
	or (I believe)
		floatplus z = get_cols(fopen("filename", "r"));

	INPUT:
	f		- input file.
	ndat		- number of elements in dat.
	dat[i]		- element to be matched to comments; see parse_file
			  and xdat below.
	dflt		- value returned in z.ddat[i] (see below) if
			  there is no match with dat[i].
	sep		- separator; see xdat below.

			  basically, if, say, the input file contains
			  the line "# noise_sheaf_length=17",
			  dat[3]="noise_sheaf_length" and sep = '=',
			  then z.ddat[3] will equal 17.  if the input
			  file contains the line "# noise_sheaf_length
			  17", again dat[3]="noise_sheaf_length", but
			  this time sep = ' ' (the default), then
			  again z.ddat[3] will equal 17.
	min_elements	- only lines with at least min_elements will
			  be saved.

	RETURNS:
	floatplus z:
	z.x[i][j]	- jth element in line i
	z.m[i]		- number of elements on line i
	z.n		- number of lines.  returns 0 if the file
			  does not exist.
	z.ddat		- data from comment; see above.  if dat[i] is
			  not found, z.ddat[i] comes back equal to dflt.
	ndat		- number of elements in z.dat.
*/
	// ---declarations
	struct floatplus z;

	// ---check that file exists.  if not, return with z.n = 0.
	//    reserve space just in case.
	if (!f)
	{
		z.n=0;
		z.m = new int[1];
		z.x = new double*[1];
		z.x[0] = new double[1];
		z.ndat=1;
		z.ddat = new double[1];
		return z;
	}

	// reserve space and make a few declarations
	ostrstream x_ptr;
	ostrstream m_ptr;

	z.ddat = parse_file(f, x_ptr, m_ptr, z.n, ndat, dat, dflt,
		sep, min_elements);

	z.m = (int*) m_ptr.str();

	// ---temporarily store the whole file
	double* x = (double*) x_ptr.str();

	// ---transfer to z.x
	z.x = newdouble(z.n, z.m);

	int cnt=0;
	for (int i=0; i < z.n; i++) for (int j=0; j < z.m[i]; j++)
		z.x[i][j] = x[cnt++];

	delete [] x;
	fclose(f);

	return z;
}

int* header_info(FILE* f, int n, char** s, int dflt=-1)
{
/*
	get info from header in ".dat" files.  assume header has the
	form

	# ...
	...
	# string=n
	...

	where string is a string and n is an integer.  stops searching
	after all strings are found, or if at the first line _not_
	beginning with #.

	INPUT:
	n	- number of strings.
	s[i]	- a list of strings, i = 0, n-1.
	dflt	- if string is not found, then dflt is returned.

	RETURN:
	info[i]	- the info pointed to by string.  e.g., if
		  s[3]="noise_random_seed" and the .dat file had the
		  line "noise_random_seed=5", then info[3] would come
		  back equal to 5.  if noise_random_seed was _not_ in
		  the .dat file, then info[3] would come back set to
		  dflt.
*/
	int* info = init_constant(n, dflt);
	char* line = new char[MLL_lib];

	while (fgets(line, MLL_lib, f))
	{
		// parse
		int n_elements;
		int* nw = parse_c(line, n_elements);

		// ---if there is at least 1 element on the line, check
		//    that first element is '#'.  if not, return.
		if (n_elements > 0) if (line[nw[0]] != '#') return info;

		// only consider lines with 2 or more elements.
		if (n_elements > 1)
		{
			// ---search for equal sign in second element.
			for (int i=0; ; i++)
			{
				if (line[nw[1]+i] == '\0') break;
				else if (line[nw[1]+i] == '=')
				{
					line[nw[1]+i] = '\0';
		
					// ---look for match
					for (int j=0; j < n; j++)
						if (!strcmp(s[j], &line[nw[1]]))
							info[j]=atoi(&line[nw[1]+i+1]);
		
					break;
				}
			}
		}
		
		delete [] nw;

		// ---see if we are all done
		int k;
		for (k=0; k < n; k++) if (info[k] == dflt) break;
		if (k == n) break;
	}

	delete [] line;

	return info;
}

double* parse_file(
	FILE* f, ostrstream& x_ptr, ostrstream& m_ptr, int& lines,
	int ndat, char** dat, double dflt, char sep=' ', int min_elements=0)
{
/*
	parse file

	INPUT:
	f		- input file.
	ndat		- number of elements in dat.
	dat[i]		- element to be matched to comments; see parse_file
			  and xdat below.
	dflt		- value returned in z.ddat[i] (see below) if
			  there is no match with dat[i].
	sep		- separator; see xdat below.

			  basically, if, say, the input file contains
			  the line "# noise_sheaf_length=17",
			  dat[3]="noise_sheaf_length" and sep = '=',
			  then z.ddat[3] will equal 17.  if the input
			  file contains the line "# noise_sheaf_length
			  17", again dat[3]="noise_sheaf_length", but
			  this time sep = ' ' (the default), then
			  again z.ddat[3] will equal 17.
	min_elements	- only lines with at least min_elements will
			  be saved.

	OUTPUT:
	x_ptr	- pointer to data.
	m_ptr	- point to number of elements on each line.

	RETURNS:
	ddat	- see comments in get_cols above.
*/
	char* line = new char[MLL_lib];
	int line_no=0, n_elements;
	lines=0;
	int cnt=0;
	int* cdat = init_constant((ndat < 1 ? 1 : ndat), 0);
	double* ddat = init_constant((ndat < 1 ? 1 : ndat), dflt);

	while (fgets(line, MLL_lib, f))
	{
		line_no++;

		// parse
		int* nw = parse_c(line, n_elements);

		// check that there is at least one element on the line
		if (n_elements == 0)
		{
			if (n_elements >= min_elements)
			{
				lines++;
				m_ptr.write((char*) &n_elements,
					sizeof(n_elements));
			}
		}
		else if (line[nw[0]] != '#')
		{
			// data exists and is not a comment
			for (int i=0; i < n_elements; i++)
			{
				double y = atof(&line[nw[i]]);
				x_ptr.write((char*) &y, sizeof(y));
			}

			lines++;
			m_ptr.write((char*) &n_elements, sizeof(n_elements));
		}
		else if (line[nw[0]] == '#' && cnt < ndat)
		{
			cnt += xdat(ndat, cdat, dat, sep, n_elements, nw, line,
				ddat);
		}

		delete [] nw;
	}

	delete [] line;

	return ddat;
}

int xdat(int ndat, int* cdat, char** dat, char sep,
	int n_line, int* nw, char* line, double* ddat)
{
/*
	see if any of the elements in line appear in dat.

	INPUT:
	ndat	- number of elements in dat.
	cdat[i]	- if 1, i^th element has been found; if 1, it has not.
		  once an element is found, it is not searched for again.
	dat[i]	- i^th elements.
	sep	- separater.  two cases:
		  sep = ' ':	compare each element in line to each
		  		element in dat until a match is found
				with dat[i].  set ddat[i] to next
				element in line.  if no next element,
				set ddat[i] to 0.
		  sep != ' ':	say sep = '='.  compare each element
		  		in line up to '=' until a match is
				found with dat[i].  set ddat[i] to the
				element after '='.  if nothing is
				there, set ddat[i] to zero.
	n_line	- number of elements in line.
	nw	- list of pointers to elements in line.
	line	- contains elements to be compared to dat.  see sep above.

	OUTPUT:
	ddat	- see sep above.

	RETURNS:
	cnt	- 1 if a match is found, 0 if no match.
*/
	int cnt=0;

	// ---loop over elments in line
	for (int i=0; i < n_line; i++)
	{
		int j, found=1;
		if (sep != ' ')
		{
			// ---parse to sep and put in '\0'.
			for (j=0;; j++)
			{
				if (line[nw[i]+j] == '\0')
				{
					found=0;
					break;
				}

				if (line[nw[i]+j] == sep)
				{
					line[nw[i]+j] = '\0';
					found=1;
					break;
				}
			}

		}

		// ---loop over lements in dat
		if (found) for (int n=0; n < ndat; n++)
		{
			if (cdat[n] == 0) if (!strcmp(dat[n], &line[nw[i]]))
			{
				cdat[n]=1;
				cnt=1;
				if (sep == ' ')
				{
					if (i+1 < n_line)
						ddat[n]=atof(&line[nw[i+1]]);
					else ddat[n]=0.0;
				}
				else
				{
					if (line[nw[i]+j+1] != '\0')
						ddat[n]=atof(&line[nw[i]+j+1]);
					else ddat[n]=0.0;
				}
			}
		}
	}
	
	return cnt;
}

int isnumber(char* s, double& x)
{
/*
	If s is a legal number, the number is transfered to x and a 1 is
	returne. Otherwise a 0 is returned.
*/
	// ---the strategy is to cycle through the elements of s and
	//    determine whether each element is legal. what is legal
	//    depends on what has come before, since numbers can
	//    includ things like 1.6e-2.
	int flag=1;
	for (int i=0; i < strlen(s); i++)
	{
		if (flag == 1)
		{
			// ---allowed values are +/-, 0-9, .
			if (s[i] == '+' || s[i] == '-')		flag=2;
			else if (isdigit(s[i]))			flag=3;
			else if (s[i] == '.')			flag=6;
			else					flag=0;
		}
		else if (flag == 2)
		{
			// ---allowed values are 0-9, .
			if (isdigit(s[i]))			flag=3;
			else if (s[i] == '.')			flag=6;
			else					flag=0;
		}
		else if (flag == 3)
		{
			// ---allowed values are e/E, 0-9, .
			if (s[i] == 'e' || s[i] == 'E')		flag=4;
			else if (isdigit(s[i]))			flag=3;
			else if (s[i] == '.')			flag=7;
			else					flag=0;
		}
		else if (flag == 4)
		{
			// ---allowed values are +/-, 0-9
			if (s[i] == '+' || s[i] == '-')		flag=5;
			else if (isdigit(s[i]))			flag=5;
			else					flag=0;
		}
		else if (flag == 5)
		{
			// ---allowed values are 0-9
			if (isdigit(s[i]))			flag=5;
			else					flag=0;
		}
		else if (flag == 6)
		{
			// ---allowed values are e/E, 0-9, .
			if (isdigit(s[i]))			flag=7;
			else					flag=0;
		}
		else if (flag == 7)
		{
			// ---allowed values are e/E, 0-9
			if (s[i] == 'e' || s[i] == 'E')		flag=4;
			else if (isdigit(s[i]))			flag=7;
			else					flag=0;
		}

		if (flag == 0) return 0;
	}

	// --- can't end on flag=2, since this indicates a + or - all alone.
	if (flag == 2) return 0;

	// ---if we've made it this far, s must be a number.
	x=atof(s);
	return 1;
}

double kl_divergence(int n, float* p, float* q)
{
/*
	computes sum_i p[i] log p[i]/q[i].  assumes p and q are
	normalized.  skips terms with q[i]=0.

	compute entropy of a probability distribution, p.
*/
	double kl=0;
	for (int i=0; i < n; i++)
		kl += p[i] > 0 && q[i] > 0 ? p[i]*log(p[i]/q[i]) : 0;

	return kl/log(2);
}

int length_string(char* s)
{
/*
	return length of \0 terminated string
*/
	int i=0;
	while (s[i++] != '\0');
	return i;
}

int max(int n, int* x)
{
/*
	return maximum of x
*/
	if (n == 0) return 0;

	int xmax = x[0];

	for (int i=1; i < n; i++) xmax = x[i] > xmax ? x[i] : xmax;

	return xmax;
}

float max(int n, float* x)
{
/*
	return maximum of x
*/
	if (n == 0) return 0;

	float xmax = x[0];

	for (int i=1; i < n; i++) xmax = x[i] > xmax ? x[i] : xmax;

	return xmax;
}

double max(int n, double* x)
{
/*
	return maximum of x
*/
	if (n == 0) return 0;

	double xmax = x[0];

	for (int i=1; i < n; i++) xmax = x[i] > xmax ? x[i] : xmax;

	return xmax;
}

int min(int n, int* x)
{
/*
	return minimum of x
*/
	if (n == 0) return 0;

	int xmin = x[0];

	for (int i=1; i < n; i++) xmin = x[i] < xmin ? x[i] : xmin;

	return xmin;
}

float min(int n, float* x)
{
/*
	return minimum of x
*/
	if (n == 0) return 0;

	float xmin = x[0];

	for (int i=1; i < n; i++) xmin = x[i] < xmin ? x[i] : xmin;

	return xmin;
}

double min(int n, double* x)
{
/*
	return minimum of x
*/
	if (n == 0) return 0;

	double xmin = x[0];

	for (int i=1; i < n; i++) xmin = x[i] < xmin ? x[i] : xmin;

	return xmin;
}

double* mmult(int n, int m, double** a, double* x)
{
/*
	return a*x.  a is n by m matrix, x is an m-dim vector
*/
	double* y = init_constant(n, (double) 0);
	for (int i=0; i < n; i++) for (int j=0; j < m; j++)
		y[i] += a[i][j]*x[j];

	return y;
}

double** mmult(int n, int m, int l, double** a, double** x)
{
/*
	return a*x.  a is n by m matrix, x is an m by l matrix
*/
	double** y = init_constant(n, l, (double) 0);
	for (int i=0; i < n; i++) for (int j=0; j < l; j++)
		for (int k=0; k < m; k++) y[i][j] += a[i][j]*x[j][k];

	return y;
}

char** cchar(int n1, int n2)
{
	char** x = (char**) calloc(n1, sizeof(char*));
	for (int n=0; n < n1; n++) x[n] = (char*) calloc(n2, sizeof(char));
	return x;
}

int** cint(int n1, int n2)
{
	int** x = (int**) calloc(n1, sizeof(int*));
	for (int n=0; n < n1; n++) x[n] = (int*) calloc(n2, sizeof(int));
	return x;
}

float** cfloat(int n1, int n2)
{
	float** x = (float**) calloc(n1, sizeof(float*));
	for (int n=0; n < n1; n++) x[n] = (float*) calloc(n2, sizeof(float));
	return x;
}

double** cdouble(int n1, int n2)
{
	double** x = (double**) calloc(n1, sizeof(double*));
	for (int n=0; n < n1; n++) x[n] = (double*) calloc(n2, sizeof(double));
	return x;
}

double** cdouble(int n1, int* n2)
{
	double** x = (double**) calloc(n1, sizeof(double*));
	for (int n=0; n < n1; n++) x[n] =
		(double*) calloc(n2[n], sizeof(double));
	return x;
}

void cfree(int n, char** x)
{
	for (int i=0; i < n; i++) free(x[i]);
	free(x);
}

void cfree(int n, int** x)
{
	for (int i=0; i < n; i++) free(x[i]);
	free(x);
}

void cfree(int n, float** x)
{
	for (int i=0; i < n; i++) free(x[i]);
	free(x);
}

void cfree(int n, double** x)
{
	for (int i=0; i < n; i++) free(x[i]);
	free(x);
}

char** newchar(int n1, int n2)
{
	char** x = new char*[n1];
	for (int n=0; n < n1; n++) x[n] = new char[n2];
	return x;
}

int** newint(int n1, int n2)
{
	int** x = new int*[n1];
	for (int n=0; n < n1; n++) x[n] = new int[n2];
	return x;
}

int** newint(int n1, int* n2)
{
	int** x = new int*[n1];
	for (int n=0; n < n1; n++) x[n] = new int[n2[n]];
	return x;
}

float** newfloat(int n1, int n2)
{
	float** x = new float*[n1];
	for (int n=0; n < n1; n++) x[n] = new float[n2];
	return x;
}

float** newfloat(int n1, int* n2)
{
	float** x = new float*[n1];
	for (int n=0; n < n1; n++) x[n] = new float[n2[n]];
	return x;
}

double** newdouble(int n1, int n2)
{
	double** x = new double*[n1];
	for (int n=0; n < n1; n++) x[n] = new double[n2];
	return x;
}

double** newdouble(int n1, int* n2)
{
	double** x = new double*[n1];
	for (int n=0; n < n1; n++) x[n] = new double[n2[n]];
	return x;
}

int* init_constant(int n, int constant)
{
/*
	overloaded function:  initialize n dimensional vector to value constant.
*/
	int* x = new int[n];
	for (int i=0; i < n; i++) x[i] = constant;
	return x;
}

int** init_constant(int n, int m, int constant)
{
/*
	overloaded function:  initialize n by m array to value constant.
*/
	int** x = newint(n, m);
	for (int i=0; i < n; i++) for (int j=0; j < m; j++) x[i][j] = constant;
	return x;
}

int** init_constant(int n, int* m, int constant)
{
/*
	overloaded function:  initialize n by m[i] array to value constant.
*/
	int** x = new int*[n];
	for (int j=0; j < n; j++)
	{
		x[j] = new int[m[j]];
		for (int i=0; i < m[j]; i++) x[j][i] = constant;
	}
	return x;
}

float* init_constant(int n, float constant)
{
/*
	overloaded function:  initialize n dimensional vector to value constant.
*/
	float* x = new float[n];
	for (int i=0; i < n; i++) x[i] = constant;
	return x;
}

float** init_constant(int n, int m, float constant)
{
/*
	overloaded function:  initialize n by m array to value constant.
*/
	float** x = newfloat(n, m);
	for (int i=0; i < n; i++) for (int j=0; j < m; j++) x[i][j] = constant;
	return x;
}

float** init_constant(int n, int* m, float constant)
{
/*
	overloaded function:  initialize n by m[i] array to value constant.
*/
	float** x = new float*[n];
	for (int j=0; j < n; j++)
	{
		x[j] = new float[m[j]];
		for (int i=0; i < m[j]; i++) x[j][i] = constant;
	}
	return x;
}

double* init_constant(int n, double constant)
{
/*
	overloaded function:  initialize n dimensional vector to value constant.
*/
	double* x = new double[n];
	for (int i=0; i < n; i++) x[i] = constant;
	return x;
}

double** init_constant(int n, int m, double constant)
{
/*
	overloaded function:  initialize n by m array to value constant.
*/
	double** x = newdouble(n, m);
	for (int i=0; i < n; i++) for (int j=0; j < m; j++) x[i][j] = constant;
	return x;
}

double** init_constant(int n, int* m, double constant)
{
/*
	overloaded function:  initialize n by m[i] array to value constant.
*/
	double** x = new double*[n];
	for (int j=0; j < n; j++)
	{
		x[j] = new double[m[j]];
		for (int i=0; i < m[j]; i++) x[j][i] = constant;
	}
	return x;
}

float mi(float** p, int imax, int jmax)
{
/*
	compute mutual information:

		M = sum_ij p_ij log(p_ij/px_i*py_j)
		px_i = sum_j p_ij
		py_j = sum_i p_ij
		

	i runs from 0 to imax - 1
	j runs from 0 to jmax - 1
*/

	// compute px and py
	float* px = new float[imax];
	for (int i=0; i < imax; i++)
	{
		px[i] = 0;
		for (int j=0; j < jmax; j++) px[i] += p[i][j];
	}
	float* py = new float[jmax];
	for (int j=0; j < jmax; j++)
	{
		py[j] = 0;
		for (int i=0; i < imax; i++) py[j] += p[i][j];
	}

	// ---mutual information
	float mi=0;
	// H(X)
	for (int i=0; i < imax; i++) mi -= px[i] > 0 ? px[i]*log(px[i]) : 0;

	// H(Y)
	for (int j=0; j < jmax; j++) mi -= py[j] > 0 ? py[j]*log(py[j]) : 0;

	// H(X,Y)
	for (int i=0; i < imax; i++) for (int j=0; j < jmax; j++)
		mi += p[i][j] > 0 ? p[i][j]*log(p[i][j]) : 0;

	mi /= log(2);

	delete [] px;
	delete [] py;

	return mi;
}

int rabs(int x)
{
	return x > 0 ? x : -x;
}

float rabs(float x)
{
	return x > 0 ? x : -x;
}

double rabs(double x)
{
	return x > 0 ? x : -x;
}

unsigned short* read_binary(FILE* f, int& n, short dum)
{
/*
	read binary file of chars.

	INPUT:
	f	- input file.
	dum	- specifies type (overloaded file).

	OUTPUT:
	n	- number of elements in file.  set to 0 if file does
		  not exist.

	RETURNS:
		- elements in file
*/
	// ---check that file exists.
	unsigned short* x;
	if (!f)
	{
		n=0;
		return x;
	}

	// ---find length
	fseek(f, 0, SEEK_END);
	int nbytes = ftell(f);
	n = nbytes/sizeof(short);
	rewind(f);

	// ---get data
	x = new unsigned short[n];
	fread(x, sizeof(short), n, f);

	// ---if dum < 0, reverse order of bytes
	if (dum < 0) for (int i=0; i < n; i++) x[i] = x[i]/256 + 256*(x[i]%256);

	// ---return
	return x;
}

int* read_binary(FILE* f, int& n, int dum)
{
/*
	read binary file of ints.

	INPUT:
	f	- input file.
	dum	- specifies type (overloaded file).

	OUTPUT:
	n	- number of elements in file.  set to 0 if file does
		  not exist.

	RETURNS:
		- elements in file
*/
	// ---check that file exists.
	int* x;
	if (!f)
	{
		n=0;
		return x;
	}

	// ---find length
	fseek(f, 0, SEEK_END);
	int nbytes = ftell(f);
	n = nbytes/sizeof((int) nbytes);
	rewind(f);

	// ---get data
	x = new int[n];
	fread(x, sizeof(x[0]), n, f);

	// ---return
	return x;
}

float* read_binary(FILE* f, int& n, float dum)
{
/*
	read binary file of floats.

	INPUT:
	f	- input file.
	dum	- specifies type (overloaded file).

	OUTPUT:
	n	- number of elements in file

	RETURNS:
		- elements in file
*/
	// ---check that file exists.
	float* x;
	if (!f)
	{
		n=0;
		return x;
	}

	// ---find length
	fseek(f, 0, SEEK_END);
	int nbytes = ftell(f);
	n = nbytes/sizeof((float) nbytes);
	rewind(f);

	// ---get data
	x = new float[n];
	fread(x, sizeof(x[0]), n, f);

	// ---return
	return x;
}

double* read_binary(FILE* f, int& n, double dum)
{
/*
	read binary file of doubles.

	INPUT:
	f	- input file.
	dum	- specifies type (overloaded file).

	OUTPUT:
	n	- number of elements in file

	RETURNS:
		- elements in file
*/
	// ---check that file exists.
	double* x;
	if (!f)
	{
		n=0;
		return x;
	}

	// ---find length
	fseek(f, 0, SEEK_END);
	int nbytes = ftell(f);
	n = nbytes/sizeof((double) nbytes);
	rewind(f);

	// ---get data
	x = new double[n];
	fread(x, sizeof(x[0]), n, f);

	// ---return
	return x;
}

int* read_binary(char* infile, int& n, int dum)
{
/*
	read binary file of ints.

	INPUT:
	infile	- file name
	dum	- specifies type (overloaded file).

	OUTPUT:
	n	- number of elements in file.  set to 0 if file does
		  not exist.

	RETURNS:
		- elements in file
*/
	// ---open file
	FILE *f = fopen(infile, "r");

	// ---check that file exists.
	int* x;
	if (!f)
	{
		n=0;
		return x;
	}

	// ---find length
	fseek(f, 0, SEEK_END);
	int nbytes = ftell(f);
	n = nbytes/sizeof((int) nbytes);
	rewind(f);

	// ---get data
	x = new int[n];
	fread(x, sizeof(x[0]), n, f);

	// ---return
	return x;
}

float* read_binary(char* infile, int& n, float dum)
{
/*
	read binary file of floats.

	INPUT:
	infile	- file name
	dum	- specifies type (overloaded file).

	OUTPUT:
	n	- number of elements in file

	RETURNS:
		- elements in file
*/
	// ---open file
	FILE *f = fopen(infile, "r");

	// ---check that file exists.
	float* x;
	if (!f)
	{
		n=0;
		return x;
	}

	// ---find length
	fseek(f, 0, SEEK_END);
	int nbytes = ftell(f);
	n = nbytes/sizeof((float) nbytes);
	rewind(f);

	// ---get data
	x = new float[n];
	fread(x, sizeof(x[0]), n, f);

	// ---return
	return x;
}

double* read_binary(char* infile, int& n, double dum)
{
/*
	read binary file of doubles.

	INPUT:
	infile	- file name
	dum	- specifies type (overloaded file).

	OUTPUT:
	n	- number of elements in file

	RETURNS:
		- elements in file
*/
	// ---open file
	FILE *f = fopen(infile, "r");

	// ---check that file exists.
	double* x;
	if (!f)
	{
		n=0;
		return x;
	}

	// ---find length
	fseek(f, 0, SEEK_END);
	int nbytes = ftell(f);
	n = nbytes/sizeof((double) nbytes);
	rewind(f);

	// ---get data
	x = new double[n];
	fread(x, sizeof(x[0]), n, f);

	// ---return
	return x;
}

double* read_n(FILE* f, int lines, int col, int& n, int silent=1)
{
/*
	read n lines of a file and return data in kth column.
	read from file f.  skips lines whose first element is #.

	NOTE:  does not check to see that file exists, or that col and
	       lines > 0!!!  do that in the calling program.

	INPUT:
	f	- file.
	lines	- number of lines to read.
	col	- column to extract.
	silent	- if on, don't list missing elements.

	OUTPUT:
	n	- number of lines actually used.

	RETURNS:
	x	- array containing n numbers.
	
*/

	// ---reserve space
	char s[MLL_lib];
	double* x = new double[lines];

	// ---main loop
	n=0;
	while (fgets(s, MLL_lib, f))
	{
		// ---get element from column "col".  ugly, but it works.
		int i=0, isav=0, n_col=0;
		while (s[i])
		{
			if (s[i] == '#') break;
			while (isspace(s[i])) i++;
			if (s[i] == '\0') break;
		
			while (!isspace(s[i]) && s[i] != '\0') i++;
			if (s[i] == '\0') break;
			s[i++] = '\0';
			if (++n_col == col) break;
			isav=i;
		}

		// ---store data
		if (n_col == col) x[n++] = atof(&s[isav]);
		if (n == lines) break;
	}

	return x;
}

plexon read_plexon(FILE* f, FILE* g=NULL, int sort=1)
{
/*
	read file assumed to be in plexon format.  if g != NULL, write
	header to file g.

	sort if sort=1.

*/
	// reserve space and make a few declarations
	char* line = new char[MLL_lib];
	ostrstream t_ptr, c_ptr;

	// ---write and write (if g != NULL) header
	while (fgets(line, MLL_lib, f))
	{
		if (!isdigit(line[0])) {if (g != NULL) fprintf(g, "%s", line);}
		else break;
	}

	// ---now get actual data.  break as soon as a non-digit or
	//    eof is encountered.  don't forget that "line" now contains
	//    the first line of data.
	char* s;
	struct plexon d;
	d.n=0;
	for (;;)
	{
		// ---parse line
		s=line;
		while (isdigit(*s)) s++;
		*s = '\0';
		s++;
		while (isspace(*s)) s++;

		// ---extract code and time
		int code = atoi(line);
		double time = atof(s);

		// ---store in osttstream
		c_ptr.write((char*) &code, sizeof(code));
		t_ptr.write((char*) &time, sizeof(time));
		d.n++;

		// ---exit at eof or non-digit (files can end with comments).
		if (!fgets(line, MLL_lib, f)) break;
		if (!isdigit(line[0])) break;
	}

	// ---save data
	d.c = (int*) c_ptr.str();
	d.t = (double*) t_ptr.str();

	// ---sort data.  use an unsophisticate method, since data
	//    should be almost sorted anyway.
	if (sort)
	{
		// ---sort data.  use an unsophisticate method, since
		//    data should be almost sorted anyway.

		int i=1;
		for (;;)
		{
			if (d.t[i] < d.t[i-1])
			{
				// ---switch and decrement i
				double tmp = d.t[i];
				d.t[i] = d.t[i-1];
				d.t[i-1] = tmp;
				int cmp = d.c[i];
				d.c[i] = d.c[i-1];
				d.c[i-1] = cmp;
				i -= i == 1 ? 0 : 1;
			}
			else i++;
			if (i == d.n) break;
		}
	}

	return d;
}

int readline(ifstream& infile, int nmax, int line, int& nchars, char* buf)
{
/*
	read a line of a file.  a couple things to note:
		1.  buf must be dimensioned to at least nmax+1.  this
		    is so a \0 can be stuck on the end.
		2.  it is important to read the whole line even if the
		    maximum number of characters has been reached.
		    otherwise, on the next call the rest of the
		    current line will be read.

	INPUT:
	infile	- input file
	nmax	- maximum number of characters allowed on line.  see
		  note (1) above.
	line	- line number; used to tell which line(s) have too
		  many characters.  if negative, no error message will
		  be written.

	OUTPUT:
	buf	- character array containing contents of line.
	nchars	- number of characters on a line.

	RETURNS:
	data	- 1 if this line has data one it.
		  0 if there are no more lines in the file.
*/
	char ch;
	int data=1, exceed=0;
	nchars=0;

	for (;;)
	{
		if (!infile.get(ch))
		{
			data = 0;
			break;
		}
		if (ch == '\n')
		{
			buf[nchars] = '\0';
			data = 1;
			break;
		}
		else if (nchars >= nmax) exceed = 1;
		else buf[nchars++] = ch;
	}
	if (exceed && line > 0)
		cerr << "too many characters on line " << line << endl;
//	cout << "readline:   " << buf << endl;
	return data;
}

int* parse_c(char* buf, int& n)
{
/*
	find all instances of whitespace delimited strings.

	INPUT:
	buf	- character array containing line to be parsed.
		  assumed to end in '\0'.

	OUTPUT:
	buf	- '\0' placed at end of every string.
	n	- number of strings in buf.

	RETURNS:
	nw	- strings are contained in buf[nw[i]], i=0, n-1.
*/
	// initialize
	n = 0;

	// skip past any blank spaces
	int i=0;
	while (isspace(buf[i])) i++;

	// variable length beast
	ostrstream i_ptr;

	// parse
	while (buf[i])
	{
		while (isspace(buf[i])) i++;
		if (buf[i] == '\0') break;

		i_ptr.write((char*) &i, sizeof(i));
		n++;

		while (!isspace(buf[i]) && buf[i] != '\0') i++;
		if (buf[i] == '\0') break;
		buf[i++] = '\0';
	}

	int* nw = (int*) i_ptr.str();

	return nw;
}

int parse_comline(int n, char** t, char* s, int inc, int dflt)
{
/*
	parse command line.  looks for match betwen t[i] and s.  two cases:
	inc > 0:  returns t[i+inc] if t[i]=s
	inc = 0:  returns 1 if t[i] = s
	in either case, if there is no match, returns dflt

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	inc	- inc > 0:  returns t[i+inc] if t[i]=s
		  inc = 0:  returns 1 if t[i] = s
		  inc < 0:  returns t[i+|inc|] if t[i]=s
		            no error if there is nothing after flag.
	dflt	- returns dflt if no match.

	RETURNS:
		- see inc and dflt above
*/
	int iflag = inc < 0 ? 0 : 1;
	inc = inc < 0 ? -inc : inc;

	int flag;
	for (flag=0; flag < n; flag++) {if (!strcmp(s, t[flag])) break;}
	if (flag == n) return dflt;

	// special case:  return 1 if this is a stand-alone flag.
	if (inc == 0) return 1;

	if (flag+inc < n)
	{
		// ---check for -flag
		int dash=0;
		for (int i=1; i <= inc; i++)
			if (t[flag+i][0] == '-' && isalpha(t[flag+i][1]))
				dash++;
		if (!dash) return atoi(t[flag+inc]);
	}

	// ---return default if inc > 1
	if (inc > 1) return dflt;

	// ---default:  error
	if (iflag) write_err("missing value after ", s, "; inc = ", inc, ".");
	else return dflt;
}

float parse_comline(int n, char** t, char* s, int inc, float dflt)
{
/*
	parse command line.  looks for match betwen t[i] and s.  two cases:
	inc > 0:  returns t[i+inc] if t[i]=s
	inc = 0:  returns 1 if t[i] = s
	in either case, if there is no match, returns dflt

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	inc	- inc > 0:  returns t[i+inc] if t[i]=s
		  inc = 0:  returns 1 if t[i] = s
	dflt	- returns dflt if no match.

	RETURNS:
		- see inc and dflt above
*/
	int flag;
	for (flag=0; flag < n; flag++) {if (!strcmp(s, t[flag])) break;}
	if (flag == n) return dflt;

	if (flag < 0) return dflt;

	// special case:  return 1 if this is a stand-alone flag.
	if (inc == 0) return 1;

	if (flag+inc < n)
	{
		// ---check for -flag
		int dash=0;
		for (int i=1; i <= inc; i++)
			if (t[flag+i][0] == '-' && isalpha(t[flag+i][1]))
				dash++;
		if (!dash) return atof(t[flag+inc]);
	}

	// ---return default if inc > 1
	if (inc > 1) return dflt;

	// ---default:  error
	cerr << "missing value after " << s << "; inc = " << inc << "." << endl;
	exit(1);
}

double parse_comline(int n, char** t, char* s, int inc, double dflt)
{
/*
	parse command line.  looks for match betwen t[i] and s.  two cases:
	inc > 0:  returns t[i+inc] if t[i]=s
	inc = 0:  returns 1 if t[i] = s
	in either case, if there is no match, returns dflt

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	inc	- inc > 0:  returns t[i+inc] if t[i]=s
		  inc = 0:  returns 1 if t[i] = s
	dflt	- returns dflt if no match.

	RETURNS:
		- see inc and dflt above
*/
	int flag;
	for (flag=0; flag < n; flag++) {if (!strcmp(s, t[flag])) break;}
	if (flag == n) return dflt;

	if (flag < 0) return dflt;

	// special case:  return 1 if this is a stand-alone flag.
	if (inc == 0) return 1;

	if (flag+inc < n)
	{
		// ---check for -flag
		int dash=0;
		for (int i=1; i <= inc; i++)
			if (t[flag+i][0] == '-' && isalpha(t[flag+i][1]))
				dash++;
		if (!dash) return atof(t[flag+inc]);
	}

	// ---return default if inc > 1
	if (inc > 1) return dflt;

	// ---default:  error
	cerr << "missing value after " << s << "; inc = " << inc << "." << endl;
	exit(1);
}

int* parse_comline(int n, char** t, char* s, int& nmatch)
{
/*
	NOTE:  this should be phased out!!!

	parse command line.  looks for match betwen t[i] and s.
	returns all elements after match until next "-" flag.
	n is number of elements.

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]

	OUTPUT:
	nmatch	- number of elements after flag.
		  if nmatch < 0, exit with error if number of elements
		  is not equal to -nmatch.

	RETURNS:
		- array containing elements.
*/
	int target = nmatch < 0 ? -nmatch : 0;

	// ---find match to flag
	int i1;
	for (i1=0; i1 < n; i1++) {if (!strcmp(s, t[i1])) break;}

	// ---no match
	if (i1 == n)
	{
		if (target) write_err(s, " needs ", target, " arguments.");
		nmatch = 0;
		int* x = new int[1];
		return x;
	}

	// ---find next occurrance of "-" followed by non-digit
	int i2;
	for (i2=i1+1; i2 < n; i2++)
	{
		if (t[i2][0] == '-') if (!isdigit(t[i2][1])) break;
	}

	nmatch = i2 - i1 - 1;
	int* x = new int[nmatch > 0 ? nmatch : 1];
	for (int i=0; i < nmatch; i++) x[i] = atoi(t[i+i1+1]);

	if (target && target != nmatch)
		write_err(s, " needs ", target, " arguments.");

	return x;
}

int* parse_comline(int n, char** t, char* s, int dum, int& nmatch)
{
/*
	a bit of ugliness to distinguish this overloaded function from
	float*, etc.:  need an int in second to last argument.

	parse command line.  looks for match betwen t[i] and s.
	returns all elements after match until next "-" flag.
	n is number of elements.

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	dum	- see below

	OUTPUT:
	nmatch	- number of elements after flag.
		  if dum < 0, exit with error if number of elements
		  is not equal to -dum.

	RETURNS:
		- array containing elements.
*/
	int target = dum < 0 ? -dum : 0;

	// ---find match to flag
	int i1;
	for (i1=0; i1 < n; i1++) {if (!strcmp(s, t[i1])) break;}

	// ---no match
	if (i1 == n)
	{
		if (target) write_err(s, " needs ", target, " arguments.");
		nmatch = 0;
		int* x = new int[1];
		return x;
	}

	// ---find next occurrance of "-" followed by non-digit
	int i2;
	for (i2=i1+1; i2 < n; i2++)
	{
		if (t[i2][0] == '-') if (!isdigit(t[i2][1])) break;
	}

	nmatch = i2 - i1 - 1;
	int* x = new int[nmatch > 0 ? nmatch : 1];
	for (int i=0; i < nmatch; i++) x[i] = atoi(t[i+i1+1]);

	if (target && target != nmatch)
		write_err(s, " needs ", target, " arguments.");

	return x;
}

int* parse_comline(int n, char** t, char* s, int dum, int& nmatch, char* dflt)
{
/*
	a bit of ugliness to distinguish this overloaded function:
	need an int in third to last argument.

	parse command line.  looks for match betwen t[i] and s.
	returns all elements after match until next "-" flag.  n is
	number of elements; i.e., the size of t.

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	dum	- < 0:  exit with error if number of elements is not
			equal to -dum.  this feature turned off if
			dflt is present.
		- > 0:  use default string (dflt) to fill in any
		        missing values, up to dum elements.
	dflt	- string containing defaults. if empty, none will be applied.

	OUTPUT:
	nmatch	- number of elements after flag.

	RETURNS:
		- array containing elements.
*/
	int ndum = (int) dum;
	int target = ndum < 0 ? -ndum : 0;

	// ---find match to flag
	int i1;
	for (i1=0; i1 < n; i1++) {if (!strcmp(s, t[i1])) break;}

	// ---no match
	if (ndum <= 0 && i1 == n)
	{
		if (target) write_err(s, " needs ", target, " arguments.");
		nmatch = 0;
		int* x = new int[1];
		return x;
	}

	// ---find next occurrance of "-" followed by non-digit and
	//    non-decimal point.
	int i2;
	for (i2=i1+1; i2 < n; i2++)
	{
		if (t[i2][0] == '-')
			if (!isdigit(t[i2][1]) && t[i2][1] != '.') break;
	}

	// ---compute length of array.
	nmatch = i2 - i1 - 1;
	int* x = new int[1 + (nmatch > ndum ? nmatch : ndum)];

	// ---set elements
	for (int i=0; i < nmatch; i++) x[i] = atoi(t[i+i1+1]);

	// ---defaults
	if (nmatch < ndum)
	{
		double* tmp = s_parse(ndum, dflt);
		for (int i=nmatch; i < ndum; i++) x[i] = (int) (0.5+tmp[i]);
		free(tmp);
	}

	// ---check that we we have hit the target number of elements,
	//    target=0 if dum > 0.
	if (target && target != nmatch)
		write_err(s, " needs ", target, " arguments.");

	return x;
}

float* parse_comline(
	int n, char** t, char* s, float dum, int& nmatch, char* dflt)
{
/*
	a bit of ugliness to distinguish this overloaded function from
	int*:  need a (float) 1 (or some such double) in third to last argument.

	parse command line.  looks for match betwen t[i] and s.
	returns all elements after match until next "-" flag.  n is
	number of elements; i.e., the size of t.

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	dum	- < 0:  exit with error if number of elements is not
			equal to -dum.  this feature turned off if
			dflt is present.
		- > 0:  use default string (dflt) to fill in any
		        missing values, up to dum elements.
	dflt	- string containing defaults. if empty, none will be applied.

	OUTPUT:
	nmatch	- number of elements after flag.

	RETURNS:
		- array containing elements.
*/
	int ndum = (int) dum;
	int target = ndum < 0 ? -ndum : 0;

	// ---find match to flag
	int i1;
	for (i1=0; i1 < n; i1++) {if (!strcmp(s, t[i1])) break;}

	// ---no match
	if (ndum <= 0 && i1 == n)
	{
		if (target) write_err(s, " needs ", target, " arguments.");
		nmatch = 0;
		float* x = new float[1];
		return x;
	}

	// ---find next occurrance of "-" followed by non-digit and
	//    non-decimal point.
	int i2;
	for (i2=i1+1; i2 < n; i2++)
	{
		if (t[i2][0] == '-')
			if (!isdigit(t[i2][1]) && t[i2][1] != '.') break;
	}

	// ---compute length of array.
	nmatch = i2 - i1 - 1;
	float* x = new float[1 + (nmatch > ndum ? nmatch : ndum)];

	// ---set elements
	for (int i=0; i < nmatch; i++) x[i] = atof(t[i+i1+1]);

	// ---defaults
	if (nmatch < ndum && strcmp(dflt, ""))
	{
		double* tmp = s_parse(ndum, dflt);
		for (int i=nmatch; i < ndum; i++) x[i] = tmp[i];
		free(tmp);
		nmatch=ndum;
	}

	// ---check that we we have hit the target number of elements,
	//    target=0 if dum > 0.
	if (target && target != nmatch)
		write_err(s, " needs ", target, " arguments.");

	return x;
}

double* parse_comline(
	int n, char** t, char* s, double dum, int& nmatch, char* dflt)
{
/*
	a bit of ugliness to distinguish this overloaded function from
	int*:  need a 1.0 (or some such double) in second to last argument.

	parse command line.  looks for match betwen t[i] and s.
	returns all elements after match until next "-" flag.
	n is number of elements.

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	dum	- < 0:  exit with error if number of elements is not
			equal to -dum.  this feature turned off if
			dflt is present.
		  > 0:  use default string (dflt) to fill in any
		  	missing values, up to dum elements.
	dflt	- string containing defaults.  if empty, none will be applied.

	OUTPUT:
	nmatch	- number of elements after flag.

	RETURNS:
		- array containing elements.
*/
	int ndum = (int) dum;
	int target = ndum < 0 ? -ndum : 0;

	// ---find match to flag
	int i1;
	for (i1=0; i1 < n; i1++) {if (!strcmp(s, t[i1])) break;}

	// ---no match
	if (i1 == n)
	{
		if (target) write_err(s, " needs ", target, " arguments.");
		nmatch = 0;
		double* x = new double[1];
		return x;
	}

	// ---find next occurrance of "-" followed by non-digit and
	//    non-decimal point.
	int i2;
	for (i2=i1+1; i2 < n; i2++)
	{
		if (t[i2][0] == '-')
			if (!isdigit(t[i2][1]) && t[i2][1] != '.') break;
	}

	// ---compute length of array.
	nmatch = i2 - i1 - 1;
	double* x = new double[nmatch > 0 ? nmatch : 1];

	// ---set elements
	for (int i=0; i < nmatch; i++) x[i] = atof(t[i+i1+1]);

	// ---defaults
	if (nmatch < ndum)
	{
		double* tmp = s_parse(ndum, dflt);
		for (int i=nmatch; i < ndum; i++) x[i] = tmp[i];
		free(tmp);
	}

	// ---check that we we have hit the target number of elements
	//    target=0 if dum > 0.
	if (target && target != nmatch)
		write_err(s, " needs ", target, " arguments.");

	return x;
}

char* parse_comline(int n, char** t, char* s, int inc, char* dflt)
{
/*
	parse command line.  looks for match betwen t[i] and s.  two cases:
	inc > 0:  returns t[i+inc] if t[i]=s
	inc = 0:  returns 1 if t[i] = s
	in either case, if there is no match, returns dflt

	INPUT:
	n	- number of elements on command line (argc)
	t	- list of strings on command line (argv)
	s	- string to be matched with t[i]
	inc	- inc > 0:  returns t[i+inc] if t[i]=s
		  inc = 0:  returns 1 if t[i] = s
	dflt	- returns dflt if no match.

	RETURNS:
		- see inc and dflt above
*/
	int flag;
	for (flag=0; flag < n; flag++) {if (!strcmp(s, t[flag])) break;}
	if (flag == n) return dflt;

	// special case:  return 1 if this is a stand-alone flag.
	if (inc == 0) return "1";

	if (flag+inc < n)
	{
		// ---check for -flag
		int dash=0;
		for (int i=1; i <= inc; i++) if (t[flag+i][0] == '-') dash++;
		if (!dash) return t[flag+inc];
	}

	// ---return default if inc > 1
	if (inc > 1) return dflt;

	// ---default:  error
	cerr << "missing value after " << s << "; inc = " << inc << "." << endl;
	exit(1);
}

double* s_parse(int n, char* v)
{
/*
	simple string parse:  turn the first n elements of v into
	doubles and store in x.  if there are fewer than n elements,
	repeat the last element that occurs.  if there are no
	elements, fill x with zeros.
*/
	int cnt=0, i=0, j=0;
	char *s, *t = (char*) calloc(1024, sizeof(char));
	double* x = (double*) calloc(n, sizeof(double));

	strcpy(t, v);
	s=t;

	for (;;)
	{
		while (*s && isspace(*s)) {s++; j++; i++;}
		while (*s && !isspace(*s)) {s++; j++;}
		if (i == j) break;

		*s='\0';
		x[cnt] = atof(&t[i]);
		cnt++;

		j++;
		s++;
		i=j;

		if (cnt == n) break;
	}

	// ---fill rest of array
	if (cnt == 0) {x[0]=0.0; cnt++;}
	if (cnt < n) for (int i=cnt; i < n; i++) x[i] = x[cnt-1];

	free(t);
	return x;
}

int unit_no(char* s)
{
/*
	get unit number.  three possibilities:
	1.  s is an integer.  return atoi(s).
	2.  s is an integer followed by a, b, c or d.  return
	    4*(n-1)+f(alpha) where n is the integer part and
	    f(a)=1, f(b)=2, f(c)=3, f(d)=4.
	3.  s is none of the above.  return 0.
*/
	int n = strlen(s);

	// ---check that first part of s is OK.
	for (int i=0; i < n-1; i++) if (!isdigit(s[i])) return 0;

	// ---if last character in s is digit, return atoi(s).
	if (isdigit(s[n-1])) return atoi(s);
	else
	{
		// ---translate last character via f(alpha) above; if
		//    last character is not a-d, return 0.
		int offset;
		if (s[n-1] == 'a') offset=1;
		else if (s[n-1] == 'b') offset=2;
		else if (s[n-1] == 'c') offset=3;
		else if (s[n-1] == 'd') offset=4;
		else return 0;

		// ---translate from plexon to unit format and return result.
		char* dum = new char[n];
		for (int i=0; i < n-1; i++) dum[i] = s[i];
		dum[n-1] = '\0';
		int cell_unit = 4*(atoi(dum)-1)+offset;

		delete [] dum;
		return cell_unit;
	}
}

void write_err(char* s)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s\n", s);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s, int n)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%d\n", s, n);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, int n, char* s2)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%d%s\n", s1, n, s2);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s, double z)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%f\n", s, z);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, char* s2)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s\n", s1, s2);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, char* s2, char* s3)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%s\n", s1, s2, s3);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, char* s2, int n, char* s3)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%d%s\n", s1, s2, n, s3);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, char* s2, int n1, char* s3, int n2)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%d%s%d\n", s1, s2, n1, s3, n2);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, char* s2, char* s3, int n, char* s4)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%s%d%s\n", s1, s2, s3, n, s4);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, char* s2, char* s3, int n1, char* s4, int n2)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%s%d%s%d\n", s1, s2, s3, n1, s4, n2);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_err(char* s1, char* s2, char* s3, char* s4, char* s5)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%s%s%s\n", s1, s2, s3, s4, s5);
	fprintf(stderr, "[00m");
	exit(1);
}

void write_warning(char* s)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s\n", s);
	fprintf(stderr, "[00m");
}

void write_warning(char* s1, char* s2, char* s3)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%s\n", s1, s2, s3);
	fprintf(stderr, "[00m");
}

void write_warning(char* s1, int n, char* s2)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%d%s\n", s1, n, s2);
	fprintf(stderr, "[00m");
}

void write_warning(char* s1, int n1, char* s2, int n2, char* s3)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%d%s%d%s\n", s1, n1, s2, n2, s3);
	fprintf(stderr, "[00m");
}

void write_warning(char* s1, char* s2, char* s3, int n, char* s4)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%s%s%d%s\n", s1, s2, s3, n, s4);
	fprintf(stderr, "[00m");
}

void write_warning(char* s1, int n1, char* s2, char* s3, int n2, char* s4)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%d%s%d%s\n", s1, n1, s2, n2, s3);
	fprintf(stderr, "[00m");
}

void write_warning(char* s1, double x, char* s2)
{
	fprintf(stderr, "[01;31m");
	fprintf(stderr, "%s%f%s\n", s1, x, s2);
	fprintf(stderr, "[00m");
}

void write_message(char* s, char* t)
{
FILE* f = fopen(s, "w");

fprintf(f, t);
fclose(f);

char sys[80];
sprintf(sys, "more %s", s);
system(sys);

sprintf(sys, "rm %s", s);
system(sys);
}

void zfree(floatplus& z)
{
	for (int i=0; i < z.n; i++) delete [] z.x[i];
	delete [] z.x;
	delete [] z.m;
}
