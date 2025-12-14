
/*=============================================================================================
    
                              EQ_SECONDO_GRADO_RESOLVERS.CPP
    
    Funzioni per la risoluzione di equazioni di secondo grado del tipo: a*x^2 + b*x + c = 0
    Il calcolo delle soluzioni avviene tramite la formula risolutiva, con attenzione
    alla stabilità numerica. 

    !! ATTENZIONE !!
    Ricordarsi di aggiungere #include "NOME_HEADER.h" nel file principale del progetto.
    Includere nell'header le seguenti librerie:
    #include <cmath>
    #include <iostream>
    #include <math.h>
    using namespace std;

    Autore: Scaffardi Alesssandro
=============================================================================================*/



//--------------------------------SECONDO_ORDINE_EQUATION_SOLVER---------------------------------
/*
Riceve in input:
 - i coefficienti a, b, c dell'equazione di secondo grado
 - un array xc di dimensione 2 in cui salvare le soluzioni trovate
*/
void Second_order_equation_solver(double a, double b, double c, double* xc){
    double delta = b*b - 4*a*c;

    if(delta < 0){
        cout << "L'equazione non ha soluzioni reali." << endl;
        xc[0] = xc[1] = NAN; // Indica che non ci sono soluzioni reali
        return;
    }

    if( b >= 0){
        xc[0] = (-2*c) / (b + sqrt(delta));
        xc[1] = (-b - sqrt(delta)) / (2*a);
    } else {
        xc[0] = (-b + sqrt(delta)) / (2*a);
        xc[1] = (-2*c) / (b - sqrt(delta));
    }
    
}





