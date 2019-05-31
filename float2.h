#ifndef FLOAT2_H
#define FLOAT2_H

/*
	Peter Latham
	July, 1997

	Creates an object of the form:

		x[i][j]	- floating point array.
		m[i]	- number of elements in x[i].
		n	- number of elements in x.

	E.G., if 

		x[0] = 7 1.3 2
		x[1] =
		x[2] = 4.2 -19
		x[3] = -3 -2 1.4 18.5

	then

		m[0] = 3
		m[1] = 0
		m[2] = 2
		m[3] = 4

		n    = 4
*/
#include <fstream>
#include <cstdlib>

class float2
{
private:

public:
	int		n;
	int*		m;
	double**	x;

	// overloaded constructor.  commented out so that class float2.h can
	// be included without worrying about float2.h, I think.
//	float2();

	// destructor
	~float2()
	{
		/*
		for some reason, code bombs when this line is included.
		for (int i=0; i < n; i++) delete [] x[i];
		*/
		delete [] x;
		delete [] m;
	}


	// initializer
	void init_float2(int n, int* m, double** x);
};

#endif
