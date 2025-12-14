#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>



using namespace std;

double forward_derivative(double (*func)(double), double x, double h);
double backward_derivative(double (*func)(double), double x, double h);
double central_derivative(double (*func)(double), double x, double h);
double fourth_order_derivative(double (*func)(double), double x, double h);
double func(double x);


