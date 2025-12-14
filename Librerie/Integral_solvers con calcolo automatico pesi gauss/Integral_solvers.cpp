
/*===================================================================================================

                                            INTEGRAL_SOLVERS.CPP

In questa libreria sono implementate alcune tecniche di integrazione numerica, in particolare, metodo
dei Rettangoli, Trapezi, Simpson e Gauss.

Per ogniuna di queste funzioni quando viene richiamata bisogna passare la funzione da integrare,
l'estremo di integrazione sinistro (a), l'estremo di integrazione destro (b) e il numero di
sotto intervalli in cui si vuole dividere l'intervallo [a,b].

E' anche implementata una funzione per il calcolo di integrali multidimensionali (2D) con metodo di Gauss,
il dominio di integrazione è considerato rettangolare.

Le funzioni da integrare devono essere implementate alla fine di questa libreria.

La libreria contiene anche una funzione per il calcolo automatico dei pesi e delle radici
del polinomio di Legendre di grado Ng, necessari per il metodo di Gauss. Con le impostazioni attuali
si possono calcolare i pesi e le radici fino a Ng = 11, per aumentare questo limite bisogna
modificare la funzione "pesi" contenuta in questa libreria, aumentando N.


!!! ATTENZIONE !!!
Bisogna includere il file nell'header ponendo:   #include "NOME_HEADER.h"
Nell'header deve anche essere presente #include<math.h>



Autore:
Scaffardi Alessandro
====================================================================================================*/
#include"header.h"


//-------------------------------------- METODO DEI RETTANGOLI --------------------------------------
/*
Riceve in input:
 - puntatore alla funzione da integrare
 - estremo di integrazione sinistro a
 - estremo di integrazione destro b
 - numero di sottointervalli N
*/
double rettangoli(double (*f)(double), double a, double b, int N){

    double h = (b-a)/double(N);
    double integral = 0;

    for(int i=0; i<N; i++){
        integral = integral + f(a+i*h)*h;
    }

    return integral;
}



//------------------------------------- METODO DEI TRAPEZI ------------------------------------------
/*
Riceve in input:
 - puntatore alla funzione da integrare
 - estremo di integrazione sinistro a
 - estremo di integrazione destro b
 - numero di sottointervalli N
*/
double trapezi(double (*f)(double), double a, double b, int N){

    double h = (b-a)/double(N);
    double integral = 0;

    for(int i=0; i<N; i++){
        integral = integral + (f(a+i*h)+f(a+(i+1)*h))*h/2;
    }

    return integral;

}



//------------------------------------- METODO DI SIMPSON -------------------------------------------
/*
Riceve in input:
 - puntatore alla funzione da integrare
 - estremo di integrazione sinistro a
 - estremo di integrazione destro b
 - numero di sottointervalli N
*/
double simpson(double (*f)(double), double a, double b, int N){

    double h = (b-a)/double(N);
    double integral = 0;

    if( N % 2 != 0){
        cout << "Nella regola di Simpson N deve essere pari!!!" <<endl;
        return NAN;
    }

    for(int i=1; i<N; i=i+2){
        integral = integral + f(a+(i-1)*h)*h/3 + f(a+i*h)*4*h/3 + f(a+(i+1)*h)*h/3;
    }


    return integral;


}



//----------------------------------------- GAUSS --------------------------------------------------
/*
Riceve in input:
 - puntatore alla funzione da integrare
 - estremo di integrazione sinistro a
 - estremo di integrazione destro b
 - numero di sottointervalli N
 - grado del polinomio di Legendre Ng
*/
double gauss(double (*f)(double), double a, double b,int N ,int Ng){

    double integral = 0, somma = 0, xc = 0, dx = 0;
    //double w[3] = {8.0/9 , 5.0/9 , 5.0/9}, x[3] = {0 , sqrt(3.0/5) , -sqrt(3.0/5)};
    double w[Ng]={};
    double x[Ng]={};
    pesi(Ng,w,x);
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






// -------------------------------- METODO DI GAUSS MULTIDIMENSIONALE -------------------------------
/*
Riceve in input:
 - puntatore alla funzione da integrare
 - estremo di integrazione sinistro xa
 - estremo di integrazione destro xb
 - estremo di integrazione sinistro ya
 - estremo di integrazione destro yb
 - numero di sottointervalli N
 - grado del polinomio di Legendre Ng
*/
double gauss2d(double (*f)(double,double), double xa, double xb,double ya, double yb,int N ,int Ng){

    double somma = 0;
    //double w[3]={8.0/9 , 5.0/9 , 5.0/9}, xg[3] = {0 , sqrt(3.0/5) , -sqrt(3.0/5)};
    //double xg[4] = { -0.8611363116,  -0.3399810436,   0.3399810436,   0.8611363116 };
    //double w[4]  = {  0.3478548451,   0.6521451549,   0.6521451549,   0.3478548451 };
    double w[Ng]={};
    double xg[Ng]={};
    pesi(Ng,w,xg);
    double dy = (yb-ya)/N , dx = (xb-xa)/N, xc = 0, yc = 0 , sumx = 0 , sumy = 0;

    for(int t = 0; t < N; t++){
    for(int j = 0; j < N; j++){

        yc = ((t+1)*dx + t*dx)*0.5 + ya;
        xc = ((j+1)*dy + j*dy)*0.5 + xa;

        sumy = 0;

        for(int i = 0; i < Ng; i++){

            sumx = 0;

        for(int k = 0; k < Ng; k++){

            sumx += w[k]*f(xc + 0.5*dx*xg[k],yc + 0.5*dy*xg[i]);

        }

            sumx *= 0.5*dx;
            sumy += w[i]*sumx;

        }

            sumy *= 0.5*dy;
            somma+=sumy;
    }
    }

    return somma;
}




//----------------------------------------- FUNZIONE GENERALE PER IL CALCOLO DELLE X E PESI PER GAUSS  -------------------------------------------
/*
Riceve in input:
 . - grado polinomio di Legendre Ng
 . - puntatori ad array in cui salvare i pesi e le radici
*/
void pesi(int Ng,double *w,double *x){

    int a = -1, b = 1, N = 100, n = Ng, Nr = 0;
    double tol = 1.0e-7;
    vector <double>root;


    Nr = general_root_solvers(func,deriv,root,a,b,tol,N,n);

    for(int i = 0; i<Nr; i++){
        x[i] = root[i];
        w[i]=(2.0/((1-x[i]*x[i])*deriv(func,x[i],n)*deriv(func,x[i],n)));
    }

    cout << "Zeri polinomio di Legendre di grado "<<n<<" e relativi pesi: "<<endl;

    for(int i = 0; i<Nr; i++){
        cout<<"Radici: "<<x[i]<<"; Pesi: "<<w[i]<<endl;
    }


}



//------------------------------------ FUNZIONI DA INTEGRARE --------------------------------------

double f(double x){

    return sin(x)*exp(-x)*cos(x/2);  

}


double f_2dim(double x,double y){

    return x*x*x*x*y*y + 2*x*x*y*y - y*x*x + 2;

}

double pigreco(double x,double y){

    if(sqrt(x*x+y*y)<=1){

        return 1;

    }
    else{

        return 0;

    }

}