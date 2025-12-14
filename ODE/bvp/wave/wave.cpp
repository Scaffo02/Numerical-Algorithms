#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;
void dYdt(double r, double *Y, double *k, double s);
void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*,double), double dt, int neq, double s);
double Residual(double s);
double Bisection(double (*func)(double), double a, double b, double tol);

int main(){

    double k = 0, k_fin = 20;
    int Nsteps = 5;
    int dk = (k_fin-k)/Nsteps; 
    vector <double>zeros; 

    for(int i = 0 ; i < Nsteps ;i++){

        double zero = Bisection(Residual,k + i*dk,k + (i+1)*dk,1.e-7);
        zeros.push_back(zero);
        cout << zero<<endl;
    }

    ofstream file("Soluzioni.dat");
    for(int i = 0; i<Nsteps;i++){

        int neq = 2;
        int Nsteps = 100;
        double r_fin = 1;
        double r = 0;
        double dr = r_fin / Nsteps;
        double s = zeros[i];
        double Y[2]={0,1};
        file << r <<"\t"<<Y[0]<<endl;
        for(int i=0; i<Nsteps; i++){
            
            RK4(r, Y, dYdt, dr, neq, s);
            r+=dr;
            file << r <<"\t"<<Y[0]<<endl;
        } 
        file <<endl<<endl;

    }
    

    file.close();



    return 0;
}



//------------------------------------------ RK4 ------------------------------------------------

void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*,double), double dt, int neq, double s){
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

    RHSFunc(t,Y,k1,s);
    for(int i = 0; i<neq ; i++){
        Y1[i] = Y[i] + dt*k1[i]*0.5;
    }

    RHSFunc(t+dt/2,Y1,k2,s);
    for(int i = 0; i<neq ; i++){
        Y2[i] = Y[i] + dt*k2[i]*0.5;
    }

    RHSFunc(t+dt/2,Y2,k3,s);
    for(int i = 0; i<neq ; i++){
        Y3[i] = Y[i] + dt*k3[i];
    }
    
    RHSFunc(t+dt,Y3,k4,s);
    for(int i = 0; i<neq ; i++){
        Y[i] = Y[i] + (dt/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    }


}



//------------------------------------ RESIDUAL -------------------------------
double Residual(double s){

    int neq = 2;
    int Nsteps = 100;
    double r_fin = 1;
    double r = 0;
    double dr = r_fin / Nsteps;

    double Y[2]={0,1};
    for(int i=0; i<Nsteps; i++){

        RK4(r, Y, dYdt, dr, neq, s);
        r+=dr;
    }  
    return Y[0] - 0;

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

void dYdt(double r, double *Y, double *k,double s){

    double x = Y[0], y = Y[1];
    k[0] = y;
    k[1] = -s*s*x;
    
}