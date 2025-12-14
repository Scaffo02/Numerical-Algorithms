#include "header.h"


int main(){

    double x = 0;
    double v = 0;
    double a = 0;
    double alfa = 10;
    int N = 1000;
    double h = alfa / N;
    differentiation(x,v,a,h,N,func,alfa);

    return 0;
}