#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


double EuleroStep(double t, double *Y, void (*RSHFunc)(double, double*, double*), double dt, int neq);
void dYdt(double t, double *Y, double *R);



