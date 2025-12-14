#include <iostream>
#include<math.h>
#include <cmath>
#include <vector>
#include <iomanip>

using  namespace std;

double rettangoli(double (*f)(double), double a, double b, int N);
double trapezi(double (*f)(double), double a, double b, int N);
double simpson(double (*f)(double), double a, double b, int N);
double gauss(double (*f)(double), double a, double b,int N ,int Ng);
void pesi(int Ng,double *w,double *x);
double gauss2d(double (*f)(double,double), double xa, double xb,double ya, double yb,int N ,int Ng);
double f(double x);
double f_2dim(double x,double y);
double pigreco(double x,double y);

double general_root_solvers(double (*func)(double,int),double (*deriv)(double (*)(double, int), double, int),vector<double>& root, double a, double b, double tol,int N,int n);
double Bisection(double (*func)(double,int), double a, double b, double tol,int n);
double FalsePos(double (*func)(double,int), double a, double b, double tol,int n);
double Secant(double (*func)(double,int), double a, double b, double tol,int n);
double Newton(double (*func)(double, int),double (*deriv)(double (*)(double, int), double, int), double tol, int n, double xc);

double func(double x, int n);
double deriv(double(*func)(double,int),double x, int n);
