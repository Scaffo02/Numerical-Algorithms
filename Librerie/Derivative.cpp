


// forward_derivative(func, t, h)
// func : funzione x(t)
// t    : punto in cui valutare la derivata
// h    : passo
double forward_derivative(double (*func)(double), double t, double h) {
    return (func(t + h) - func(t)) / h;
}



// backward_derivative(func, t, h)
// func : funzione x(t)
// t    : punto in cui valutare la derivata
// h    : passo
double backward_derivative(double (*func)(double), double t, double h) {
    return (func(t) - func(t - h)) / h;
}



// central_derivative(func, t, h)
// func : funzione x(t)
// t    : punto in cui valutare la derivata
// h    : passo
double central_derivative(double (*func)(double), double t, double h) {
    return (func(t + h) - func(t - h)) / (2 * h);
}



// fourth_order_derivative(func, t, h)
// Derivata centrale con accuratezza O(h^4)
// func : funzione x(t)
// t    : punto in cui valutare la derivata
// h    : passo
double fourth_order_derivative(double (*func)(double), double t, double h) {
    return ( -func(t + 2*h) + 8*func(t + h) - 8*func(t - h) + func(t - 2*h) ) 
           / (12 * h);
}



// second_derivative(func, t, h)
// Seconda derivata numerica (schema centrale)
// func : funzione x(t)
// t    : punto in cui valutare la seconda derivata
// h    : passo
double second_derivative(double (*func)(double), double t, double h) {
    return (func(t + h) - 2*func(t) + func(t - h)) / (h * h);
}



// func(t)
// t : variabile reale
// Ritorna un esempio di funzione a scelta (qui puoi mettere quella che vuoi)
double func(double t) {
    return sin(t);  // esempio; puoi modificarla con quella che ti serve
}
