#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;

void JacobiMethod(double **psi1, double **psi2, double **S, double h, int NX, int NY);
void GaussSiedel(double **psi, double **S, double h, int NX, int NY);
void SOR(double **psi, double **S, double h, int NX, int NY, double omega);


int main() {

    int NX = 32, NY = 32;
    double h = 1.0 / (NX - 1);
    double tol = 1e-7;
    double omega = 2.0/(1+M_PI/NX);

    int k_J = 0;
    int k_GS = 0;
    int k_SOR = 0;

    // Allocate matrices
    double **psi1 = new double*[NX];
    double **psi2 = new double*[NX];
    double **psi = new double*[NX];
    double **psi_sor = new double*[NX];
    double **S = new double*[NX];

    psi1[0] = new double[NX * NY];
    psi2[0] = new double[NX * NY];
    psi[0] = new double[NX * NY];
    psi_sor[0] = new double[NX * NY];
    S[0] = new double[NX * NY];

    for(int i=1;i<NX;i++){
        psi[i]  = psi[i-1] + NY;
        psi1[i] = psi1[i-1] + NY;
        psi2[i] = psi2[i-1] + NY;
        psi_sor[i] = psi_sor[i-1] + NY;
        S[i]    = S[i-1] + NY;
    }

    // Set S(x,y)=0 (Laplace equation, change if needed)
    for(int i=0;i<NX;i++)
        for(int j=0;j<NY;j++)
            S[i][j] = 0.0;

    // Apply boundary conditions
    for(int i=0;i<NX;i++){
        for(int j=0;j<NY;j++){

            double x = i*h;
            double y = j*h;

            // Bordo: funzione analitica f(x,y)
            if(i==0 || i==NX-1 || j==0 || j==NY-1){
                psi1[i][j] = exp(-M_PI*x) * sin(M_PI*y)+S[i][j]*(x*x + y*y)/4.0;
                psi2[i][j] = psi1[i][j];
                psi[i][j]  = psi1[i][j];
                psi_sor[i][j] = psi1[i][j];
            }
            else{
                psi_sor[i][j] = 0.0;
                psi[i][j]  = 0.0;
                psi1[i][j] = 0.0;
                psi2[i][j] = 0.0;
            }
        }
    }

    // Jacobi iterations
    double err;
    do{
        err = 0.0;

        JacobiMethod(psi1, psi2, S, h, NX, NY);

        for(int i=1;i<NX-1;i++){
            for(int j=1;j<NY-1;j++){
                err += fabs(psi2[i+1][j]-2*psi2[i][j]+psi2[i-1][j]+psi2[i][j+1]-
                            2*psi2[i][j]+psi2[i][j-1]-h*h*S[i][j]);
                psi1[i][j] = psi2[i][j];
            }
        }
        k_J ++;

    } while(err > tol);

    // Gauss-Seidel iteration
    do{
        err = 0.0;

        GaussSiedel(psi, S, h, NX, NY);

        for(int i=1;i<NX-1;i++){
            for(int j=1;j<NY-1;j++){
                err += fabs(psi[i+1][j]-2*psi[i][j]+psi[i-1][j]+psi[i][j+1]-
                            2*psi[i][j]+psi[i][j-1]-h*h*S[i][j]);
                }
        }

        k_GS ++;
    } while(err > tol);


    // SOR iteration
    do{
        err = 0.0;

        SOR(psi_sor, S, h, NX, NY, omega);

        for(int i=1;i<NX-1;i++){
            for(int j=1;j<NY-1;j++){
                err += fabs(psi_sor[i+1][j]-2*psi_sor[i][j]+psi_sor[i-1][j]+psi_sor[i][j+1]-
                            2*psi_sor[i][j]+psi_sor[i][j-1]-h*h*S[i][j]);
                }
        }

        k_SOR ++;
    } while(err > tol);


    // Save results to file
    ofstream fout("solution.dat");
    for(int i=0;i<NX;i++){
        for(int j=0;j<NY;j++){
            double x = i*h;
            double y = j*h;
            fout << x << " " << y << " " << psi1[i][j] <<" "<< psi[i][j] <<" "<<psi_sor[i][j]<< "\n";
        }
        fout << "\n";
    }
    fout.close();

    cout <<"Step per raggiungere la convergenza:"<<endl;
    cout <<"Jacobi: "<< k_J <<endl;
    cout <<"Gauss-Siedel: "<< k_GS <<endl;
    cout <<"SOR: "<< k_SOR <<endl;


    // Free memory
    delete[] psi_sor[0]; delete[] psi_sor;
    delete[] psi[0]; delete[] psi;
    delete[] psi1[0]; delete[] psi1;
    delete[] psi2[0]; delete[] psi2;
    delete[] S[0]; delete[] S;

    return 0;
}


void JacobiMethod(double **psi1, double **psi2, double **S, double h, int NX, int NY){

    for(int i=1;i<NX-1;i++){
        for(int j=1;j<NY-1;j++){
            psi2[i][j] = 0.25 * (psi1[i+1][j] + psi1[i-1][j]
                               + psi1[i][j+1] + psi1[i][j-1]
                               - h*h*S[i][j]);
        }
    }
}


void GaussSiedel(double **psi, double **S, double h, int NX, int NY){

    for(int i=1;i<NX-1;i++){
        for(int j=1;j<NY-1;j++){
            psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-h*h*S[i][j]);            
        }
    }

}

void SOR(double **psi, double **S, double h, int NX, int NY, double omega){
    
    for(int i=1;i<NX-1;i++){
        for(int j=1;j<NY-1;j++){
            psi[i][j] = (1-omega)*psi[i][j]+(omega/4)*(psi[i-1][j]+psi[i+1][j]+psi[i][j-1]+psi[i][j+1]-h*h*S[i][j]);            
        }
    }
}