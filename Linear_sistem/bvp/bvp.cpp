#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;
void solve_tridiagonal_full(double  *a, double *b, double *c,double *r,double *x, int n);

int main(){

    ofstream file("Soluzione.dat");
    int n = 52;
    double x_f = 1.0, x0 = 0.0;
    double dx = (x_f-x0)/double(n-1);
    double a[n] = {};
    double b[n] = {};
    double c[n] = {};
    double r[n] = {};
    double x[n] = {};

    for(int i = 0; i<n ; i++){
        a[i] = 1.0;
        b[i] = -2.0;
        c[i] = 1.0;
        if(i == 1){
            r[i] = 1.0*dx*dx - 0.0;
        }else if(i == n-2){
            r[i] = 1.0*dx*dx - 0.0;
        }else{
            r[i] = 1.0*dx*dx;
        }

    }
    
    solve_tridiagonal_full(a+1,b+1,c+1,r+1,x+1,n-2);
    double s = 0;
    for(int i = 0; i<n; i++){
        s = i*dx +x0; 
        file << s << "\t" << x[i] <<endl;
        
    }

    file.close();
   


    return 0;
}

void solve_tridiagonal_full(double  *a, double *b, double *c,double *r,double *x, int n){

    double h[n-1] = {};
    double p[n]   = {};


    for(int i = 0; i<n-1 ; i++){
        if(i == 0){
            h[i] = c[i]/b[i];
        }else{
            h[i] = c[i]/(b[i]-a[i]*h[i-1]);
        }
    }

    for(int i = 0; i<n ; i++){
        if(i == 0){
            p[i] = r[i]/b[i];
        }else{
            p[i] = (r[i]-a[i]*p[i-1])/(b[i]-a[i]*h[i-1]);
        }
    }
   
    x[n-1] = p[n-1];
    for(int i = n-2; i>=0 ; i--){
        x[i] = p[i] - h[i]*x[i+1];
    }

}
