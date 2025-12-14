
/*==================================================================================================

                                        ODE_RESOLVERS.CPP   



Autore:
Scaffardi Alessandro

====================================================================================================*/




#include "header.h"


//------------------------------------ GENERAL RESOLVERS --------------------------------------------

void ODE_resolvers(int neq, double *Y, double dt, double t, int Nsteps){
    ofstream file("Soluzione.dat");
    int i = 0;

    cout <<"Scegli il metodo da usare: \n 1) Eulero \n 2) RK_2 \n 3) RK_4\n";
    cin >> i;


    switch(i){

        case 1:
        
        for (int i = 0; i < Nsteps; i++) {
            file << t << "\t" << Y[0] << "\t" << Y[1] << endl;
            EuleroStep(t, Y, dYdt, dt, neq);
            t += dt;
        }
        break;


        case 2 : 

        for (int i = 0; i < Nsteps; i++) {
            file << t << "\t" << Y[0] << "\t" << Y[1] << endl;
            RK2(t, Y, dYdt , dt, neq);
            t += dt;
        }
        break;

        case 3:

        for (int i = 0; i < Nsteps; i++) {
            file << t << "\t" << Y[0] << "\t" << Y[1] << endl;
            RK4(t, Y, dYdt , dt, neq);
            t += dt;
        }
        break;

        default:
        cout <<"Immissione non valida!!"<<endl;

    }

        

    file.close();
}






//------------------------------------ METODO DI EULERO ----------------------------------------------

void EuleroStep(double t, double *Y, void (*RSHFunc)(double, double*, double*), double dt, int neq){
    // Take one step dt using Euler method for the solution of dY/dt = rhs. 
    // Here neq is the number of ODE (the dimensionality of Y[]) and *RHSFunc() is a 
    // pointer to the actual function (in this case it points to dYdt()) that 
    // calculates the right-hand side of the system of equations.
    int k = 0;
    double rhs[256];

    RSHFunc(t, Y, rhs);
    for(k=0; k<neq; k++){
        Y[k] += dt*rhs[k];
    }
}



//------------------------------------------ RK2 ------------------------------------------------

void RK2(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){

    double Y1[neq], k1[neq], k2[neq];

    RHSFunc(t,Y,k1);
    for(int i = 0; i<neq ; i++){
        Y1[i] = Y[i] + 0.5*dt*k1[i]; 
    }

    RHSFunc(t+0.5*dt,Y1,k2);
    for(int i = 0; i<neq ; i++){
        Y[i] += dt*k2[i];
    }

}



//------------------------------------------ RK4 ------------------------------------------------

void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){

    double Y1[neq], Y2[neq], Y3[neq], k1[neq], k2[neq], k3[neq], k4[neq];

    RHSFunc(t,Y,k1);
    for(int i = 0; i<neq ; i++){
        Y1[i] = Y[i] + dt*k1[i]*0.5;
    }

    RHSFunc(t+dt/2,Y1,k2);
    for(int i = 0; i<neq ; i++){
        Y2[i] = Y[i] + dt*k2[i]*0.5;
    }

    RHSFunc(t+dt/2,Y2,k3);
    for(int i = 0; i<neq ; i++){
        Y3[i] = Y[i] + dt*k3[i];
    }
    
    RHSFunc(t+dt,Y3,k4);
    for(int i = 0; i<neq ; i++){
        Y[i] = Y[i] + (dt/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    }


}



//-------------------------------------- ODE -----------------------------------------------

void dYdt(double t, double *Y, double *k){

    double x = Y[0], y = Y[1];
    k[0] = y;
    k[1] = -x;
    
}



