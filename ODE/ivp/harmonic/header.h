#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

void ODE_resolvers(int neq, double *Y, double dt, double t, int Nsteps);
void RK2(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
void velocity_velet(double *x, double *v, double h, int neq);
void position_velet(double *x, double *v, double dt,int neq);
void ODE_resolvers(int neq, double t, double dt, double *Y,int Nsteps);
void EuleroStep(double t, double *Y, void (*RSHFunc)(double, double*, double*), double dt, int neq);
void acc(double *x,double *a,int neq);
void dYdt(double t, double *Y, double *k);

const double omega = 1;
