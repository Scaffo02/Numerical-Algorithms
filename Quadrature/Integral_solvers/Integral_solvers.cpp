#include "header_integral.h"



double rettangoli(double (*f)(double), double a, double b, int N){

    double h = (b-a)/double(N);
    double integral = 0;

    for(int i=0; i<N; i++){
        integral = integral + f(a+i*h)*h;
    }

    return integral;
}

double trapezi(double (*f)(double), double a, double b, int N){

    double h = (b-a)/double(N);
    double integral = 0;

    for(int i=0; i<N; i++){
        integral = integral + (f(a+i*h)+f(a+(i+1)*h))*h/2;
    }

    return integral;

}

double simpson(double (*f)(double), double a, double b, int N){

    double h = (b-a)/double(N);
    double integral = 0;

    for(int i=1; i<N; i=i+2){
        integral = integral + f(a+(i-1)*h)*h/3 + f(a+i*h)*4*h/3 + f(a+(i+1)*h)*h/3;
    }


    return integral;


}

double gauss(double (*f)(double), double a, double b,int N ,int Ng){

    double integral = 0, somma = 0, xc = 0, dx = 0;
    double w[3] = {8.0/9 , 5.0/9 , 5.0/9}, x[3] = {0 , sqrt(3.0/5) , -sqrt(3.0/5)};
    double h = (b-a)/N;


    for(int t = 0; t < N; t++){

        somma = 0;
        dx = (((t+1)*h+a)-(t*h+a))/2;
        xc = (((t+1)*h+a)+(t*h+a))/2;

        for(int i = 0; i < Ng; i++){

            somma = somma + w[i]*f(dx*x[i] + xc);

        }

        integral = integral + dx*somma;

    }

    return integral;
}


double f(double x){

    return x*x-sin(x);

}
