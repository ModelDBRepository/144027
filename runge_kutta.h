#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

//  headers
#include <stdlib.h>
#include "derivs.h"
#include "float2.h"

//  definitions

//  prototypes
class runge_kutta
{
private:
	float **dx1, **dx2, **dx3, **dx4;
	float2 xint;

public:
	// constructor
	runge_kutta(int n, int* m);
	runge_kutta(int n, int mm);

	// fourth order Runge-Kutta
	void rk4(derivs& dxdt, float2& x, float& t, float dt);

	// second order Runge-Kutta
	void rk2(derivs& dxdt, float2& x, float& t, float dt);
};

#endif

