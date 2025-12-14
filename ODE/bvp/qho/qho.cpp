#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;
void dYdt(double x, double *Y, double *k);
void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq);
double Residual(double xm);
double Bisection(double (*func)(double), double a, double b, double tol);
double general_root_solvers(double (*func)(double),vector<double>& root, double a, double b, double tol,int N);


int main(){

    vector<double>root;
    double Ea = 0;
    double Eb = 5;
    double tol = 1.0e-10;
    int N = 50;

    double Ne = general_root_solvers(Residual,root,Ea,Eb,tol,N);
   

    return 0;
}







//-------------------------------------- RSIDUAL -----------------------------------------------
 
double Residual(double E){

    double xm = 0.3;

    // Calcoliamo il forward
    double x = -10.0;
    int Nsteps = 10000;
    double dx = (xm-x)/Nsteps;
    int neq = 3;
    double Y[neq]={exp(-x*x/2),-exp(-x*x/2)*x,E};
    
    for(int i = 0; i < Nsteps; i++){

        RK4(x,Y,dYdt,dx,neq);
        x = x + dx;
        
    }
    double yL = Y[0];
    double dyL = Y[1];
    
    //Calcoliamo backword
    x = 10.0;
    dx = (xm-x)/Nsteps;

    Y[0]=exp(-x*x/2);
    Y[1]=-exp(-x*x/2)*x;
    for(int i = 0; i < Nsteps; i++){

        RK4(x,Y,dYdt,dx,neq);
        x = x + dx;
        
    }
    double yR = Y[0];
    double dyR = Y[1];

    //Calcoliamo il residuo
    double res = (dyL*yR - dyR*yL)/(sqrt(dyL*yR*dyL*yR+dyR*yL*dyR*yL) + 1.0e-10);
    return res;

}



//------------------------------------------ BRACKETING ---------------------------------------------------------
/*
Riceve in input:
 - puntatori a funzione func e deriv
 - vettore root in cui salvare le radici trovate
 - estremi dell'intervallo [a,b]
 - tolleranza tol
 - numero di suddivisioni N dell'intervallo [a,b]
*/
double general_root_solvers(double (*func)(double),vector<double>& root, double a, double b, double tol,int N){

    int Nr = 0;
    int nr = 0;
    int i = 0;
    double dx = (b-a)/N;
    vector<double>xL;
    vector<double>xR;
    vector<double>root_Bisezione;

    double x_L = 0;
    double x_R = 0;

    for(int i = 0; i < N-1 ; i++){
        x_L = a + i*dx;
        x_R = a + (i+1)*dx;

        if(func(x_L)*func(x_R)<=0){
            xL.push_back(x_L);
            xR.push_back(x_R);
            Nr++;
        }

    }


    for(int i = 0; i < Nr; i++){

        double controllo = Bisection(func,xL[i],xR[i],tol);
        if(i==0){
            root_Bisezione.push_back(controllo);
            cout <<"Radice: "<<root_Bisezione.back()<<endl;
            nr++;

        }
        if(controllo != root_Bisezione.back() && i!=0){
            root_Bisezione.push_back(controllo);
            cout <<"Radice: "<<root_Bisezione.back()<<endl;
            nr++;
        }

    }



 
    for(int i = 0; i < nr; i++)
    {
        root.push_back(root_Bisezione[i]);
    }
      


    return nr;


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




//-------------------------------------- ODE ----------------------------------

void dYdt(double x, double *Y, double *k){

    double psi = Y[0], y = Y[1], E = Y[2];
    k[0] = y;
    k[1] = -2*(E*psi-0.5*x*x*psi);
    
}