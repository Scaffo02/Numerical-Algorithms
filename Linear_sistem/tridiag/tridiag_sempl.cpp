#include <iostream>
using namespace std;

void solve_tridiagonal_full(double  *a, double *b, double *c,double *r,double *x, int n);


int main(){

    int n = 5;
    double a[n] = {1,1,1,1,1}; //il primo valore non viene utilizzato
    double b[n]   = {2,2,2,2,2};
    double c[n] = {1,1,1,1,1}; //l'ultimo valore non viene utilizzato
    double r[n]   = {1,0,3,1,0};
    double x[n]   = {0};

    solve_tridiagonal_full(a,b,c,r,x,n);
    cout << "soluzione sistema A:" <<endl;
    for(int i = 0; i<n; i++){
        cout << x[i] <<endl;
    }

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