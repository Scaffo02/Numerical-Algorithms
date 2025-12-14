

#include "header.h"




int main(){


    double integral = 0;
    int n = 0;
    cout <<setprecision(10);


    cout <<"Inserisci il grado del polinomio di Legendre da utilizzare per integrare:"<<endl;
    cin >>n;

    integral = gauss(f,0,10,100,n);

    cout <<"Valore integrale 1 dim:"<<integral<<endl;

    cout <<"==========================================================="<<endl;
    cout <<"==========================================================="<<endl;


   integral = gauss2d(pigreco,-1,1,-1,1,1000,n);
   cout <<"Valore integrale 2 dim:"<<integral<<endl;



    return 0;
}
