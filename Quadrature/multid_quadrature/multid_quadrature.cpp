#include <iostream>
#include <cmath>
#include <iomanip>



double f(double x,double y);
double f2(double x,double y);
double gauss2d(double (*f)(double,double), double xa, double xb,double ya, double yb,int N ,int Ng);
using namespace std;



int main()
{
    cout << setprecision(5);
    cout <<setiosflags(ios:: scientific);

    int Ng = 4, N = 10;
    double integ_gauss_pol = 0, integ_gauss_cil = 0, integ_gauss_da_1D_pol = 0, xa = -1, xb = 1, ya = -1, yb = 1;

    cout <<"Calcolo funzione plinomiale 2d con metodo di Gauss tra -1<x<1 e -1<y<1:"<<endl;

    integ_gauss_pol = gauss2d(f,xa,xb,ya,yb,N,Ng);

    cout <<"Gauss 2D genuino: "<<integ_gauss_pol<<endl;
    cout <<"Esatto: " << 412.0/45<<endl;


    cout <<"Calcolo pigreco di Gauss:"<<endl;
    do{
        integ_gauss_cil = gauss2d(f2,xa,xb,ya,yb,N,Ng);
        N = N*2;
    }while(abs(integ_gauss_cil- M_PI)>pow(10,-5));
    cout <<"Gauss (N = "<<N<<"): "<<integ_gauss_cil<<endl;
    cout <<"Esatto: " << M_PI<<endl;

    return 0;
}





// Questa funzione calcola direttamente l'integrale in 2D con il metodo di Gauss con Ng = 4
double gauss2d(double (*f)(double,double), double xa, double xb,double ya, double yb,int N ,int Ng){

    double somma = 0;
    //double w[3]={8.0/9 , 5.0/9 , 5.0/9}, xg[3] = {0 , sqrt(3.0/5) , -sqrt(3.0/5)};
    double xg[4] = { -0.8611363116,  -0.3399810436,   0.3399810436,   0.8611363116 };
    double w[4]  = {  0.3478548451,   0.6521451549,   0.6521451549,   0.3478548451 };
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









double f(double x,double y){

    return x*x*x*x*y*y + 2*x*x*y*y - y*x*x + 2;

}


double f2(double x,double y){

    if(sqrt(x*x+y*y)<=1){

        return 1;

    }
    else{

        return 0;

    }

}




