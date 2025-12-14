#include "header.h"


double forward_derivative(double (*func)(double), double x, double h) {
    return (func(x + h) - func(x)) / h;
}

double backward_derivative(double (*func)(double), double x, double h) {
    return (func(x) - func(x - h)) / h;
}

double central_derivative(double (*func)(double), double x, double h) {
    return (func(x + h) - func(x - h)) / (2 * h);
}

double fourth_order_derivative(double (*func)(double), double x, double h) {
    return ( -func(x + 2*h) + 8*func(x + h) - 8*func(x - h) + func(x - 2*h) ) / (12 * h);
}   

double func(double x) {
    return sin(x); // Example function: f(x) = sin(x)
}