


/*==============================================================================================================

                                                ROOT_SOLVERS.CPP

Questa librearia contiene alcune funzioni per il calcolo degli zeri di una funzione. In particolare sono implementati
i metodi di BISEZIONE, FALSE POSITION, SECANTE e NEWTON. In fondo alla libreria bisogner� definire la funzione e la sua
derivata per cui si vuole calcolare gli zeri su un dato intervallo.
Ricordare che questi metodi prevedo che nell'intervallo considerato ci sia al pi� uno zero.

Se si vuole calcolare tutti gli zeri su un dato intervallo per una funzione si pu� richiamare la funzione
GENERAL_ROOT_SOLVERS che permette di trovare in generale gli zeri di una funzione su un determinato intervallo [a,b].
Il programma in se non trova gli zeri ma suddivide l'intervallo in N intervalli su cui si utilizzaranno le funzioni
descritte in precedenza, i risultati vengono salvati e sono quindi disponibili per essere analizzati
successivamente.

Gli algoritmi hanno una soglia di iterazione di 500 passi di default che si pu� modificare.


!!! ATTENZIONE !!!
Bisogna includere il file nell'header ponendo:

#include "NOME_HEADER.h"



Autore:
Scaffardi Alessandro
================================================================================================================*/




//------------------------------------------ BRACKETING ---------------------------------------------------------
/*
Riceve in input:
 - puntatori a funzione func e deriv
 - vettore root in cui salvare le radici trovate
 - estremi dell'intervallo [a,b]
 - tolleranza tol
 - numero di suddivisioni N dell'intervallo [a,b]
*/
double general_root_solvers(double (*func)(double),double (*deriv)(double),vector<double>& root, double a, double b, double tol,int N){

    int Nr = 0;
    int nr = 0;
    int i = 0;
    double dx = (b-a)/N;
    vector<double>xL;
    vector<double>xR;
    vector<double>root_Bisezione;
    vector<double>root_FalsPos;
    vector<double>root_Secante;
    vector<double>root_Newton;

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


    cout <<"Metodo della Bisezione:" <<endl;
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

    cout <<"----------------------------------------------------------"<<endl;

    cout <<"Metodo dei Falsi Positivi:" <<endl;
    for(int i = 0; i < Nr; i++){
        double controllo = FalsePos(func,xL[i],xR[i],tol);
        if(i == 0){
            root_FalsPos.push_back(controllo);
            cout <<"Radice: "<<root_FalsPos.back()<<endl;

        }
        if(controllo != root_FalsPos.back() && i!=0){
            root_FalsPos.push_back(controllo);
            cout <<"Radice: "<<root_FalsPos.back()<<endl;

        }
    }

    cout <<"----------------------------------------------------------"<<endl;

    cout <<"Metodo della secante:" <<endl;
    for(int i = 0; i < Nr; i++){
        double controllo = Secant(func,xL[i],xR[i],tol);
        if( i== 0 ){
            root_Secante.push_back(controllo);
            cout <<"Radice: "<<root_Secante.back()<<endl;

        }
        if(controllo != root_Secante.back() && i!=0){
            root_Secante.push_back(controllo);
            cout <<"Radice: "<<root_Secante.back()<<endl;

        }
    }

    cout <<"----------------------------------------------------------"<<endl;

    cout <<"Metodo di Newton:" <<endl;
    for(int i = 0; i < Nr; i++){
        double xc = (xL[i]+xR[i])/2;
        double controllo = Newton(func,deriv,xc,tol);
        if( i == 0 ){
            root_Newton.push_back(controllo);
            cout <<"Radice: "<<root_Newton.back()<<endl;

        }
        if(controllo != root_Newton.back() && i!=0){
            root_Newton.push_back(controllo);
            cout <<"Radice: "<<root_Newton.back()<<endl;

        }
    }


    cout <<"----------------------------------------------------------"<<endl;
    cout <<" Di quale metodo vuoi salvare il risultato?"<<endl;
    cout <<" 1) Bisezione"<<endl;
    cout <<" 2) False Position"<<endl;
    cout <<" 3) Secante"<<endl;
    cout <<" 4) Newton"<<endl;
    cin >>i;


    switch(i){

        case 1:
            for(int i = 0; i < nr; i++)
            {
                root.push_back(root_Bisezione[i]);
            }
        break;

        case 2:
            for(int i = 0; i < nr; i++)
            {
                root.push_back(root_FalsPos[i]);
            }
        break;

        case 3:
            for(int i = 0; i < nr; i++)
            {
                root.push_back(root_Secante[i]);
            }
        break;

        case 4:
            for(int i = 0; i < nr; i++)
            {
                root.push_back(root_Newton[i]);
            }
        break;

        default:
            cout <<"Valore inserito non ammissibile!!"<<endl;
    }


    return nr;


}



//---------------------------------- METODO DELLA BISEZIONE -----------------------------------------------------
/*
Riceve in input:
 - puntatore a funzione func
 - estremi dell'intervallo [a,b]
 - tolleranza tol
*/
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



//---------------------------------- METODO FALSE POSITION  -----------------------------------------------------
/*
Riceve in input:
 - puntatore a funzione func
 - estremi dell'intervallo [a,b]
 - tolleranza tol
*/
double FalsePos(double (*func)(double), double a, double b, double tol){

    double xm = -func(a)*((b-a)/(func(b)-func(a)))+a;
    double F_xm = func(xm);
    double xm_0 = xm;
    double del = a-xm;
    int k = 1;

    if(func(a) == 0){
        return a;
    }
    if(func(b) == 0){
        return b;
    }

    while( fabs(del) > tol ){

        if(F_xm == 0){
            return xm;
        }

        if(func(a)*F_xm<0){

            b = xm;

        }

        if(func(a)*F_xm>0){

            a = xm;

        }

        xm = -func(a)*((b-a)/(func(b)-func(a)))+a;
        F_xm = func(xm);
        del = xm-xm_0;
        xm_0 = xm;
        k = k+1;


        if(k>500){
            cout <<"!!Numero massimo di iterazioni dell'algoritmo raggiunto (default 500), se si vuole continuare modificare k nella funzione contenuta nel file root_solvers.cpp!!"<<endl;
            break;
        }

    }

    return xm;

}



//---------------------------------- METODO DELLA SECANTE -----------------------------------------------------
/*
Riceve in input:
 - puntatore a funzione func
 - estremi dell'intervallo [a,b]
 - tolleranza tol
*/
double Secant(double (*func)(double), double a, double b, double tol){

    double Fa = func(a), Fb = func(b);
    double dx = b-a;
    int k = 1;

    if(func(a) == 0){
        return a;
    }
    if(func(b) == 0){
        return b;
    }

    while(fabs(dx)>tol){

        dx = Fb*(b-a)/(Fb-Fa);

        a = b;
        Fa = Fb;
        b = b-dx;
        Fb = func(b);
        k++;

        if(k>500){
            cout <<"!!Numero massimo di iterazioni dell'algoritmo raggiunto (default 500), se si vuole continuare modificare k nella funzione contenuta nel file root_solvers.cpp!!"<<endl;
            break;
        }

        if(fabs(dx)>100){
            cout <<"!!L'algoritmo sta divergendo!!"<<endl;
            return nan("");
        }

    }

    return b;

}



//---------------------------------- METODO DI NEWTON ----------------------------------------------------------
/*
Riceve in input:
 - puntatori a funzione func e deriv
 - xc punto iniziale
 - tolleranza tol
*/
double Newton(double (*func)(double), double (*deriv)(double), double xc, double tol){

    double xc_1 = 0;
    double dx = xc_1 - xc;
    int k = 1;

    if(func(xc) == 0){
        return xc;
    }


    do{

        xc_1 = xc - func(xc)/deriv(xc);

        dx = xc_1 - xc;
        xc = xc_1;
        k++;

        if(k>500){
            cout <<"!!Numero massimo di iterazioni dell'algoritmo raggiunto (default 500), se si vuole continuare modificare k nella funzione contenuta nel file root_solvers.cpp!!"<<endl;
            break;
        }

    }while(fabs(dx) > tol);

    return xc;

}




//---------------------------------- FUNZIONE E RISPETTIVA DERIVATA SU CUI VALUTARE GLI ZERI  -----------------------------------------------------

double func(double x){

    return sin(x);

}

// La derivata serve per il metodo di Newton
double deriv(double x){

    return cos(x) ;

}



