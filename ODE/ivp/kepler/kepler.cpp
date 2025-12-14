#include "header.h"


int main(){

    /*Il numero di equazioni è da impostare pari a:
        - 4 se si usano i metodi RK o Eulero (1,2,3)
        - 2 se si usano i metodi di Velet (3,4)
    */
    int neq = 2;  

    // Condizioni iniziali dimensionali (cgs)
    double x0_cgs = 1.5e13;   // cm, 1 AU
    double y0_cgs = 0.0;
    double r0_cgs = sqrt(x0_cgs*x0_cgs + y0_cgs*y0_cgs);
    double v0_cgs = sqrt(GM/r0_cgs); // velocità orbitale [cm/s]

    double vx0_cgs = 0;
    double vy0_cgs = v0_cgs;

    // Conversione in unità adimensionali
    double alfa = 1;
    double x0 = x0_cgs / L_ref;
    double y0 = y0_cgs / L_ref;
    double vx0 = vx0_cgs * T_ref / L_ref*sqrt(alfa);
    double vy0 = vy0_cgs * T_ref / L_ref*sqrt(alfa);

    double Y[4] = {x0, y0, vx0, vy0};

    double dt = 0.01;           // passo adimensionale
    double t_in = 0.0;
    double t_fin = 100.0;        // tempo adimensionale
    int Nsteps = int((t_fin-t_in)/dt);

    ODE_resolvers(neq, Y, dt, t_in, Nsteps);

    return 0;
}