#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void ODE_resolvers(int neq, double *Y, double dt, double t, int Nsteps);
void EuleroStep(double t, double *Y, void (*RSHFunc)(double, double*, double*), double dt, int neq);
void RK2(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
void dYdt(double t, double *Y, double *R);
