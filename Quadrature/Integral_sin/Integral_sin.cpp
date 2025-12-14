#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>


double simpson(double (*f)(double), double a, double b, int N);
double trapezi(double (*f)(double), double a, double b, int N);
double gauss(double (*f)(double), double a, double b, int N,int Ng);
double f(double x);
double f2(double x);


using namespace std;


int main()
{
    cout << setprecision(6);
    cout <<setiosflags(ios:: scientific);


    int Ng = 3, N = 1;
    double integ_simp = 0 , integ_gauss = 0 , integ_trap = 0, a = 0, b=0;


    cout <<"Calcolo sin(x)/x da [0,x], inserisci x:"<<endl;
    cin >> b;


    do{

        integ_trap = trapezi(f,a,b,N);
        cout <<"Trapezi con n = "<< N <<":"<<integ_trap<<endl;

        N = N*2;

    }while(N<=8);


    cout <<"-----------------------------------------------------"<<endl;


    N = 2;

    do{

        integ_simp = simpson(f,a,b,N);
        cout <<"Simpson con n = "<<N<<":"<<integ_simp<<endl;

        N = N*2;

    }while(N<=8);


    cout <<"-----------------------------------------------------"<<endl;


    N = 1;

    integ_gauss = gauss(f,a,b,N,Ng);
    cout <<"Gauss con n = "<<N<<":"<<integ_gauss<<endl;


    double n = 10;
    double x_i = 1.0/n;
    double x[int((b-a)*n)] = {0}, integ[int((b-a)*n)] = {0};
    N = (b-a)*2;


    for(int i = 0; i < n*(b-a); i++ ){

        x[i] = a + x_i*i;
        integ[i] = gauss(f,a,x[i],N,Ng);

    }


    ofstream fdata;
    fdata.open("int_sin.dat");


     for(int i = 0; i < n*(b-a); i++ ){

        fdata << x[i]<<" "<<integ[i] <<endl;

    }

    fdata.close();


    return 0;
}



double simpson(double (*f)(double), double a, double b, int N){


    double h = double(b-a)/double(N);
    double integral = 0;


    for(int i=1; i<N; i=i+2){

        integral = integral + f(a+(i-1)*h)*h/3 + f(a+i*h)*4*h/3 + f(a+(i+1)*h)*h/3;

    }


    return integral;
}

double trapezi(double (*f)(double), double a, double b, int N){


    double h = double(b-a)/double(N);
    double integral = 0;


    for(int i=0; i<N; i++){

        integral = integral + (f(a+i*h)+f(a+(i+1)*h))*h/2;

    }


    return integral;
}



double gauss(double (*f)(double), double a, double b,int N ,int Ng){

    double integral = 0, somma = 0, xc = 0, dx = 0;
    double w[3] = {8.0/9 , 5.0/9 , 5.0/9}, x[3] = {0 , sqrt(3.0/5) , -sqrt(3.0/5)};
    double h = double(b-a)/N;


    for(int t = 0; t < N; t++){

        somma = 0;
        dx = (((t+1)*h+a)-(t*h+a))/2;
        xc = (((t+1)*h+a)+(t*h+a))/2;

        for(int i = 0; i < Ng; i++){

            somma = somma + w[i]*f(dx*x[i] + xc);

        }

        integral = integral + dx*somma;

    }

    return integral;
}


double f(double x){

    if (x != 0){

        return sin(x)/x;

    }
    else{

        return 1;

    }
}


