// ======================================================
//      Derivate numeriche di una traiettoria x(t)
// ======================================================
// funzioni:
//  - differentiation: calcola e salva x(t), v(t), a(t)
//  - forward_derivative: derivata in avanti
//  - backward_derivative: derivata all’indietro
//  - central_derivative: derivata centrale
//  - fourth_order_derivative: derivata centrale 4º ordine
//  - second_derivative: seconda derivata
//  - func: funzione x(t) da derivare

#include "header.h"

// differentiation(x, v, a, h, N, func, alfa)
// x, v, a : variabili posizione, velocità, accelerazione (inizialmente ignorate)
// h       : passo temporale
// N       : numero di passi (interi, ma passato come double)
// func    : puntatore alla funzione x(t, alfa)
// alfa    : parametro della funzione x(t)
// Output: file "trajectory.dat" con t, x, v, a
void differentiation(double x, double v , double a, double h, double N,
                     double (*func)(double,double), double alfa){
    ofstream outfile("trajectory.dat");

    for (int i = 0; i <= N; i++) {
        double t = i * h;
        if(t == 0){
            x = func(t, alfa);                              // posizione x(t)
            double v_0 = forward_derivative(func, t, h, alfa);       // velocità dx/dt
            v   = forward_derivative(func, t+h,h ,alfa);    
            a = (v-v_0)/h;        // accelerazione d²x/dt²
            outfile << fixed << setprecision(6) << t << " " 
                << x << " " << v_0 << " " << a << endl;
        }
        else{
            x = func(t, alfa);                              // posizione x(t)
            v = central_derivative(func, t, h, alfa);       // velocità dx/dt
            a = second_derivative(func, t, h, alfa);        // accelerazione d²x/dt²
            outfile << fixed << setprecision(6) << t << " " 
                << x << " " << v << " " << a << endl;
            x = 0.0;
            v = 0.0;
            a = 0.0;
        }
        
    }

   outfile.close();
}

// forward_derivative(func, t, h, alfa)
// func : funzione x(t, alfa)
// t    : punto in cui valutare la derivata
// h    : passo
// alfa : parametro di func
double forward_derivative(double (*func)(double,double), double t, double h, 
                          double alfa) {
    return (func(t + h, alfa) - func(t, alfa)) / h;
}

// backward_derivative(func, t, h, alfa)
// func : funzione x(t, alfa)
// t    : punto in cui valutare la derivata
// h    : passo
// alfa : parametro di func
double backward_derivative(double (*func)(double,double), double t, double h, 
                           double alfa) {
    return (func(t, alfa) - func(t - h, alfa)) / h;
}

// central_derivative(func, t, h, alfa)
// func : funzione x(t, alfa)
// t    : punto in cui valutare la derivata
// h    : passo
// alfa : parametro di func
double central_derivative(double (*func)(double,double), double t, double h, 
                          double alfa) {
    return (func(t + h, alfa) - func(t - h, alfa)) / (2 * h);
}

// fourth_order_derivative(func, t, h, alfa)
// Derivata centrale con accuratezza O(h^4)
// func : funzione x(t, alfa)
// t    : punto in cui valutare la derivata
// h    : passo
// alfa : parametro di func
double fourth_order_derivative(double (*func)(double,double), double t, double h,   
                               double alfa) {
    return ( -func(t + 2*h, alfa) + 8*func(t + h, alfa) - 8*func(t - h, alfa) + 
              func(t - 2*h, alfa) ) / (12 * h);
}   

// second_derivative(func, t, h, alfa)
// Seconda derivata numerica (schema centrale)
// func : funzione x(t, alfa)
// t    : punto in cui valutare la seconda derivata
// h    : passo
// alfa : parametro di func
double second_derivative(double (*func)(double, double), double t, double h, 
                         double alfa) {
    return (func(t + h, alfa) - 2*func(t, alfa) + func(t - h, alfa)) / (h * h);
}

// func(t, alfa)
// t    : tempo
// alfa : parametro reale
// Ritorna x(t) = alfa t^2 − t^3 (1 − exp(−alfa^2 / t))
double func(double t,double alfa) {
  
    if( t == 0){
        return 0;
    }
    return alfa * t * t - t * t * t * (1 - exp(-alfa * alfa / t)); 
}
