
#include "header.h"


int main(){
    ofstream file;
    double h = 0.5;
    double x = 1.0;
    double err_central = 0.0, err_forward = 0.0, err_backward = 0.0, err_fourth_order = 0.0;
    double exact = cos(x);

    file.open("errors.txt");
    for (int i = 0; i < 10; i++) {
        err_central = fabs(central_derivative(func, x, h) - exact);
        err_forward = fabs(forward_derivative(func, x, h) - exact);
        err_backward = fabs(backward_derivative(func, x, h) - exact);
        err_fourth_order = fabs(fourth_order_derivative(func, x, h) - exact);
        h /= 2;
        file << 1.0/h << " " << err_forward << " " << err_backward << " " << err_central << " " << err_fourth_order << endl;
    }

    file.close();


}
   

