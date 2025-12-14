#include "header.h"



int main(){

    ofstream file("Traiettoria.txt");
    int neq = 1;          // Number of equations
    double Y[neq] = {1.0};  // Initial conditions
    double dt = 0.01;      // Time step
    double t = 0.0;      // Initial time
    int Nsteps = 500;     // Number of steps


    for (int i = 0; i < Nsteps; i++) {
        printf("t = %.2f, y = %.6f\n", t, Y[0]);
        file << t << "\t" << Y[0] << endl;
        EuleroStep(t, Y, dYdt, dt, neq);
        t += dt;
    }

    return 0;
}