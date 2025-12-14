#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
void dYdt(double r, double *Y, double *k);
void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
double Residual(double s);
double Bisection(double (*func)(double), double a, double b, double tol);

int main(){

    ofstream file("Soluzione.dat");

    double a = 0, b = 5, tol = 1.e-7;
    double s = Bisection(Residual,a,b,tol);
    cout << s <<endl;



    int neq = 2;
    int Nsteps = 1000;
    double r_fin = 20;
    double r = 0;
    double dr = r_fin / Nsteps;
    double Y[2]={0,s};

    file << r << "\t"<<Y[1]<<endl;

    for(int i=0; i<Nsteps; i++){

        RK4(r, Y, dYdt, dr, neq);
        r+=dr;
        file << r << "\t"<<Y[0]<<endl;
    }

    file.close();

    return 0;
}



//------------------------------------------ RK4 ------------------------------------------------

void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){
/**
  Esegue un singolo passo di integrazione con il metodo di Runge-Kutta del 4° ordine (RK2)
  per il sistema dY/dt = RHS.
 
  - t         Tempo corrente.
  - Y         Vettore delle variabili del sistema (modificato in uscita).
  - RHSFunc   Puntatore alla funzione che calcola il termine destro del sistema (es. dYdt()).
  - dt        Passo temporale di integrazione.
  - neq       Numero di equazioni del sistema (dimensione del vettore Y[]).
 */
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



//------------------------------------ RESIDUAL -------------------------------
double Residual(double s){

    int neq = 2;
    int Nsteps = 1000;
    double r_fin = 20;
    double r = 0;
    double dr = r_fin / Nsteps;

    

    
    double Y[2]={0,s};
    r = 0;
    for(int i=0; i<Nsteps; i++){

        RK4(r, Y, dYdt, dr, neq);
        r+=dr;

        }
        
    return Y[0] - 1;

}



//-------------------------- BISEZIONE ----------------------------------------
double Bisection(double (*func)(double), double a, double b, double tol){

    double dx = b-a;
    double xm = (a+b)/2;
    double F_xm = func(xm);
    int k = 1;

    if(func(a) == 0){
        return a;
    }
    if(func(b) == 0){
        return b;
    }

    while( fabs(dx) > tol ){

        if(F_xm == 0){

            return xm;

        }

        if(func(a)*F_xm<0){

            b = xm;

        }

        if(func(a)*F_xm>0){

            a = xm;

        }

        dx = b-a;
        xm = (a+b)/2;
        F_xm = func(xm);
        k = k+1;

        if(k>500){
            cout <<"!!Numero massimo di iterazioni dell'algoritmo raggiunto (default 500), se si vuole continuare modificare k nella funzione contenuta nel file root_solvers.cpp!!"<<endl;
            break;
        }

    }

    return xm;

}

//-------------------------------------- ODE ----------------------------------

void dYdt(double r, double *Y, double *k){

    double x = Y[0], y = Y[1];
    double rho = exp(-r)/(8*M_PI);
    k[0] = y;
    k[1] = -4*M_PI*r*rho;
    
}