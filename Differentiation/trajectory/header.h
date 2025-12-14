#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

double forward_derivative(double (*func)(double,double), double x, double h, double alfa);
double backward_derivative(double (*func)(double,double), double x, double h, double alfa);
double central_derivative(double (*func)(double,double), double x, double h, double alfa);
double fourth_order_derivative(double (*func)(double,double), double x, double h, double alfa);
double second_derivative(double (*func)(double,double), double x, double h, double alfa);
double func(double x,double alfa);

void differentiation(double x, double v , double a, double h, double N, double (*func)(double,double), double alfa);

