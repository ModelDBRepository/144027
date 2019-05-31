/*

	Peter Latham
	July, 1997

*/

#include "float2.h"

// constructor
void float2::init_float2(int n, int* m, double** x)
{
	int i, j;

	// set up class variables
	this->m = new int[n];
	this->x = new double*[n];
	for (i=0; i < n; i++) this->x[i] = new double[m[i]];

	// transfer data
	this->n = n;
	for (i=0; i < n; i++)
	{
		this->m[i] = m[i];
		for (j=0; j < m[i]; j++) this->x[i][j] = x[i][j];
	}
}
