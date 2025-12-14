#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void ODE_resolvers(int neq, double *Y, double dt, double t, int Nsteps);
void EuleroStep(double t, double *Y, void (*RSHFunc)(double, double*, double*), double dt, int neq);
void RK2(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
void position_velet(double *x, double *v, double dt,int neq);
void velocity_velet(double *x, double *v, double dt,int neq);
void dYdt(double t, double *Y, double *R);
void acc(double *x,double *a,int neq);


//------------------- COSTANTI DI RIFERIMENTO (CGS) ------------------------
const double L_ref = 1.5e13;      // cm, esempio: 1 AU
const double GM = 1.327e26;       // cm^3/s^2, G*M Sole
const double T_ref = 3.15e7;      // s, esempio: 1 anno
//const double T_ref = sqrt((L_ref*L_ref*L_ref)/GM);    
