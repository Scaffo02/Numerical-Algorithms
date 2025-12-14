#include "header.h"

int main(){

    /*
    Il numero di equazioni deve essere impostato in base al metodo utilizzato:
    - 2 se si usano i metodi 1,2,3
    - 1 se si usano i metodi 4,5
    */
    int neq = 1;

    double T = double(2*M_PI)/omega;
    double h = 0.005*T;
    int Nsteps = 1000;
    double t = 0;
    double Y[2] = {1,0};

    ODE_resolvers(neq,t,h,Y,Nsteps);

    return 0;
}


