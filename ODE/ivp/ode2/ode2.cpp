#include "header.h"


int main(){

    int neq = 2;
    double Y[neq] = {1.0, 0.0};
    double dt = 0.01;
    double t_fin = 20*M_PI, t_in = 0;
    int Nsteps = int((t_fin-t_in)/dt);


    ODE_resolvers( neq , Y , dt , t_in , Nsteps);


    return 0;
}