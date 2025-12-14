
/*=========================================================================================================

                                        MONTECARLO_INTEGRAL.CPP

Questo file contiene l'implementazione della funzione monteCarloIntegral, che calcola l'integrale di una funzione
utilizzando il metodo Monte Carlo. La funzione accetta come parametri un puntatore a funzione, i limiti
d'integrazione x1 e x2, i limiti y1 e y2, e il numero di punti N da utilizzare per il calcolo.
La funzione restituisce il valore approssimato dell'integrale.

!! ATTENZIONE !!
Bisogna includere il file nell'header ponendo:
#include "NOME_HEADER.h"

Nell'header bisogna inserire le librerie:
#include <cstdlib>
#include <fstream>


Autore:
Scaffardi Alessandro
=========================================================================================================*/

#include "header.h"



//----------------------------------monteCarloIntegral---------------------------------------------------
/*
Riceve in input:
 - un puntatore a funzione f(double)
 - i limiti di integrazione x1 e x2
 - i limiti y1 e y2
 - il numero di punti N da utilizzare per il calcolo
*/
double monteCarloIntegral(double (*f)(double), double x1, double x2,double y1, double y2, int N) {
    double I;
    double A = x2 - x1;
    double B = y2 - y1;
    double area = A*B;
    double x, y;
    int N_in = 0;

    srand48(time(NULL));
    ofstream file("monte_carlo_results.dat");
    for(int i = 0; i < N; i++){
        x = drand48()*(x2 - x1) + x1;
        y = drand48()*(y2 - y1) + y1;

            if(y >= 0 && y <= f(x) && f(x) >= 0) {
                N_in = N_in +1;
                file << x << " " << y << endl;
            }
            else if (y < 0 && y >= f(x) && f(x) < 0) {
                N_in = N_in - 1;
                file << x << " " << y << endl;
            }
    }

    I = area*double(N_in)/N;
    file.close();
    return I;
}

double f(double x) {
    return sin(x);
}