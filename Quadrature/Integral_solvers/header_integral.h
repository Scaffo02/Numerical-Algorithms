#include<iostream>
#include<cmath>
#include<iomanip>

using namespace std;

double rettangoli(double (*f)(double), double a, double b, int N);
double trapezi(double (*f)(double), double a, double b, int N);
double simpson(double (*f)(double), double a, double b, int N);
double gauss(double (*f)(double), double a, double b,int N ,int Ng);
double f(double x);
