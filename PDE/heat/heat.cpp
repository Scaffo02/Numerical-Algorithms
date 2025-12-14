#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm> 
using namespace std;

// Funzioni
void Jacobi(double **psi, double **psi_new, double **S, double h, int NX, int NY);
double GaussSeidel(double **psi, double **S, double h, int NX, int NY);
double SOR(double **psi, double **S, double h, int NX, int NY, double omega);
void ApplyNeumann(double **psi, double hx, int NX, int NY);

int main() {
    int NX = 129, NY = 65;
    double x0=0, xf=2, y0=0, yf=1;
    double hx = (xf-x0)/(NX-1);
    double hy = (yf-y0)/(NY-1); 
    double tol = 1e-7;
    double omega = 2.0/(1 + sin(M_PI/NX)); 
    

    int kJ=0, kGS=0, kSOR=0;

    // Allocazione griglie
    double **psi_j     = new double*[NX];
    double **psi_jnew  = new double*[NX];
    double **psi_gs    = new double*[NX];
    double **psi_sor   = new double*[NX];
    double **S         = new double*[NX];

    psi_j[0]     = new double[NX*NY];
    psi_jnew[0]  = new double[NX*NY];
    psi_gs[0]    = new double[NX*NY];
    psi_sor[0]   = new double[NX*NY];
    S[0]         = new double[NX*NY];

    for(int i=1;i<NX;i++){
        psi_j[i]     = psi_j[i-1]+NY;
        psi_jnew[i]  = psi_jnew[i-1]+NY;
        psi_gs[i]    = psi_gs[i-1]+NY;
        psi_sor[i]   = psi_sor[i-1]+NY;
        S[i]         = S[i-1]+NY;
    }

    // Inizializzazione sorgente
    for(int i=0;i<NX;i++)
        for(int j=0;j<NY;j++)
            S[i][j] = 0.0;

    // Condizioni al contorno e guess iniziale
    for(int i=0;i<NX;i++){
        double x = i*hx;
        for(int j=0;j<NY;j++){
            if(j==0){      // y=0 Dirichlet
                psi_j[i][j]=psi_jnew[i][j]=psi_gs[i][j]=psi_sor[i][j]=0;
            } else if(j==NY-1){ // y=1 Dirichlet
                psi_j[i][j]=psi_jnew[i][j]=psi_gs[i][j]=psi_sor[i][j]=2-x;
            } else{        // interno
                psi_j[i][j]=psi_jnew[i][j]=psi_gs[i][j]=psi_sor[i][j]=0;
            }
        }
    }

    //---------------- Jacobi ----------------
    double err;
    do{
        err = 0;
    
        // 1. Calcola i nuovi valori per gli interni
        Jacobi(psi_j, psi_jnew, S, hx, NX, NY);
    
        // 2. Applica le condizioni di Neumann sui bordi x
        ApplyNeumann(psi_jnew, hx, NX, NY); 

        // 3. Calcola l'errore
        for(int i=1;i<NX-1;i++)
            for(int j=1;j<NY-1;j++)
                err = max(err, fabs(psi_jnew[i][j]-psi_j[i][j])); // Standard Max Norm

        // 4. Update
        for(int i=0;i<NX;i++)
            for(int j=1;j<NY-1;j++)
                psi_j[i][j] = psi_jnew[i][j];

        kJ++;
    } while(err>tol);

    cout <<"Jacobi fatto! Iter: " << kJ << endl;

    //---------------- Gauss-Seidel ----------------
    do{
        err = GaussSeidel(psi_gs,S,hx,NX,NY);
        ApplyNeumann(psi_gs, hx, NX, NY);
        kGS++;
    } while(err>tol);
    cout <<"Gauss-Seidel fatto! Iter: " << kGS << endl;

    //---------------- SOR ----------------
    do{
        err = SOR(psi_sor,S,hx,NX,NY,omega);
        ApplyNeumann(psi_sor,hx,NX,NY);
        kSOR++;
    } while(err>tol);
    cout <<"SOR fatto! Iter: " << kSOR << endl;

    
    //---------------- Save results ----------------
    ofstream fout("solution.dat");
    for(int i=0;i<NX;i++){
        double x=i*hx;
        for(int j=0;j<NY;j++){
            double y=j*hy;
            fout << x << " " << y << " " << psi_j[i][j] << " " 
                 << psi_gs[i][j] << " " << psi_sor[i][j] << "\n";
        }
        fout << "\n";
    }
    fout.close();

    // Free memory 
    delete[] psi_j[0]; delete[] psi_jnew[0]; delete[] psi_gs[0]; delete[] psi_sor[0]; delete[] S[0];
    delete[] psi_j; delete[] psi_jnew; delete[] psi_gs; delete[] psi_sor; delete[] S;

    return 0;
}

//------------------------------------
void Jacobi(double **psi, double **psi_new, double **S, double h, int NX, int NY){
    for(int i=1;i<NX-1;i++)
        for(int j=1;j<NY-1;j++)
            psi_new[i][j] = 0.25*( psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1]-h*h*S[i][j] );
}


//-----------------------------------
double GaussSeidel(double **psi, double **S, double h, int NX, int NY){
    double old_val;
    double err = 0.0; 
    for(int i=1;i<NX-1;i++){
        for(int j=1;j<NY-1;j++){
            old_val = psi[i][j];
            psi[i][j] = 0.25*(psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1] - h*h*S[i][j]); 
            err = max(err, fabs(psi[i][j]-old_val));
        }
    }
    return err;
}

//----------------------------------
double SOR(double **psi, double **S, double h, int NX, int NY, double omega){
    double old_val;
    double err = 0.0; 
    for(int i=1;i<NX-1;i++){
        for(int j=1;j<NY-1;j++){
            old_val = psi[i][j];
            double gs_part = (psi[i+1][j]+psi[i-1][j]+psi[i][j+1]+psi[i][j-1] - h*h*S[i][j]);
            
            psi[i][j] = (1.0-omega)*old_val + (omega/4.0)*gs_part;
            
            err = max(err, fabs(psi[i][j]-old_val));
        }
    }

    return err;
}


//------------------------------------
void ApplyNeumann(double **psi, double hx, int NX, int NY){
    for(int j=1;j<NY-1;j++){
        psi[0][j]     = psi[1][j];           // dψ/dx=0 
        psi[NX-1][j]  = psi[NX-2][j] + 3*hx; // dψ/dx=3 
    }
}