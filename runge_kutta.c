/*


	Peter Latham
	June, 1997

*/

#include "runge_kutta.h"

// overloaded constructor.  use the second form if m[i] is the same for all i.
runge_kutta::runge_kutta(int n, int* m)
{
	// float** variables
	this->dx1 = new float*[n];
	this->dx2 = new float*[n];
	this->dx3 = new float*[n];
	this->dx4 = new float*[n];

	// float2 variables
	this->xint.n = n;
	this->xint.m = new int[n];
	this->xint.x = new double*[n];
	for (int i=0; i < n; i++)
	{
		this->dx1[i] = new float[m[i]];
		this->dx2[i] = new float[m[i]];
		this->dx3[i] = new float[m[i]];
		this->dx4[i] = new float[m[i]];
		this->xint.x[i] = new double[m[i]];
	}

}

runge_kutta::runge_kutta(int n, int mm)
{
	// float** variables
	this->dx1 = new float*[n];
	this->dx2 = new float*[n];
	this->dx3 = new float*[n];
	this->dx4 = new float*[n];

	// float2 variables
	this->xint.n = n;
	this->xint.m = new int[n];
	this->xint.x = new double*[n];
	for (int i=0; i < n; i++)
	{
		this->dx1[i] = new float[mm];
		this->dx2[i] = new float[mm];
		this->dx3[i] = new float[mm];
		this->dx4[i] = new float[mm];
		this->xint.x[i] = new double[mm];
	}
}

// integrator
void runge_kutta::rk4(derivs& dxdt, float2& x, float& t, float dt)
{
/*
	fourth order Runge-Kutta
	INPUT:
	func   = the external subroutine used to compute values of the 
	         time derivative(s) call is func(dx,x,t)
	x      = variable to be integrated
	t      = time, the independent variable
	dt     = time step
*/
	// ---first round
	dxdt.func(dx1, x, t);
	for (int n=0; n < x.n; n++)
		for (int m=0; m < x.m[n]; m++)
			xint.x[n][m]=x.x[n][m]+0.5*dt*dx1[n][m];
	t=t+0.5*dt;

	// ---second round
	dxdt.func(dx2, xint, t);
	for (int n=0; n < x.n; n++)
		for (int m=0; m < x.m[n]; m++)
			xint.x[n][m]=x.x[n][m]+0.5*dt*dx2[n][m];

	// ---third round
	dxdt.func(dx3, xint, t);
	for (int n=0; n < x.n; n++)
		for (int m=0; m < x.m[n]; m++)
			xint.x[n][m]=x.x[n][m]+dt*dx3[n][m];
	t=t+0.5*dt;

	// ---fourth round
	dxdt.func(dx4, xint, t);
	for (int n=0; n < x.n; n++)
		for (int m=0; m < x.m[n]; m++)
			x.x[n][m]=x.x[n][m]+dt*
				(dx1[n][m]+2*(dx2[n][m]+dx3[n][m])+dx4[n][m])/6;
}

void runge_kutta::rk2(derivs& dxdt, float2& x, float& t, float dt)
{
/*
	second order Runge-Kutta
	INPUT:
	func   = the external subroutine used to compute values of the 
	         time derivative(s) call is func(dx,x,t)
	x      = variable to be integrated
	t      = time, the independent variable
	dt     = time step
*/
	int m, n;

	// dx1 = dx/dt at t, x
	dxdt.func(dx1, x, t);

	// xint.x = x+dx1*dt/2
	for (n=0; n < x.n; n++)
		for (m=0; m < x.m[n]; m++)
			xint.x[n][m]=x.x[n][m]+0.5*dt*dx1[n][m];

	// increment t by dt/2
	t=t+0.5*dt;

	// dx2 = dx/dt at t+dt/2, x+dx1*dt/2
	dxdt.func(dx2, xint, t);

	// x = x + dt*dx2
	for (n=0; n < x.n; n++)
		for (m=0; m < x.m[n]; m++) x.x[n][m]=x.x[n][m]+dt*dx2[n][m];

	// increment t by dt/2
	t=t+0.5*dt;
}
