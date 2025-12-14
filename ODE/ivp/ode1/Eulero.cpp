
#include "header.h"


double EuleroStep(double t, double *Y, void (*RSHFunc)(double, double*, double*), double dt, int neq){
    // Take one step dt using Euler method for the solution of dY/dt = rhs. 
    // Here neq is the number of ODE (the dimensionality of Y[]) and *RHSFunc() is a 
    // pointer to the actual function (in this case it points to dYdt()) that 
    // calculates the right-hand side of the system of equations.
    int k = 0;
    double rhs[256];

    RSHFunc(t, Y, rhs);
    for(k=0; k<neq; k++){
        Y[k] += dt*rhs[k];
    }

    return 0;
}

void dYdt(double t, double *Y, double *R){
    // Compute the right-hand side of the ODE dy/dt = -t*y
    double y = Y[0];
    R[0] = -t*y;
}





