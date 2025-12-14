#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
void printMatrix(double **M, int nrow, int ncol);
void eq_resolver(double **A, double * b, double *x, int n);
void eq_resolver_pp(double **A, double *b, double *x, int n);


int main(){

    int n = 3;
    int nb = 4;
    int nc = 4;
    double **A = new double*[n];
    double **B = new double*[nb];
    double **C = new double*[nc];
    double ba[n]={1,2,1};
    double bb[nb]={5,16,22,15};
    double bc[nb]={5,16,22,15};
    double xa[n];
    double xb[n];
    double xc[n];
    A[0] = new double[n*n];
    B[0] = new double[nb*nb];
    C[0] = new double[nc*nc];
    for(int i = 1; i < n; i++) A[i] = A[i-1] + n;
    for(int i = 1; i < nb; i++) B[i] = B[i-1] + nb;
    for(int i = 1; i < nb; i++) C[i] = C[i-1] + nc;


    A[0][0] = 2,  A[0][1] =-1,  A[0][2] = 0; 
    A[1][0] =-1,  A[1][1] = 2,  A[1][2] =-1; 
    A[2][0] = 0,  A[2][1] =-1,  A[2][2] = 2; 

    B[0][0] = 1,  B[0][1] = 2, B[0][2] = 1,  B[0][3] =-1; 
    B[1][0] = 3,  B[1][1] = 2, B[1][2] = 4,  B[1][3] = 4; 
    B[2][0] = 4,  B[2][1] = 4, B[2][2] = 3,  B[2][3] = 4;
    B[3][0] = 2,  B[3][1] = 0, B[3][2] = 1,  B[3][3] = 5;

    C[0][0] = 1,  C[0][1] = 2, C[0][2] = 1,  C[0][3] =-1; 
    C[1][0] = 3,  C[1][1] = 6, C[1][2] = 4,  C[1][3] = 4; 
    C[2][0] = 4,  C[2][1] = 4, C[2][2] = 3,  C[2][3] = 4;
    C[3][0] = 2,  C[3][1] = 0, C[3][2] = 1,  C[3][3] = 5;

    cout << "matrice A:" <<endl;
    printMatrix(A,n,n);
    cout << "matrice B:" <<endl;
    printMatrix(B,nb,nb);
    cout << "matrice C:" <<endl;
    printMatrix(C,nc,nc);
    eq_resolver(A,ba,xa,n);
    eq_resolver(B,bb,xb,nb);
    eq_resolver_pp(C,bc,xc,nc);
    
    cout << "soluzione sistema A:" <<endl;
    for(int i = 0; i<n; i++){
        cout << xa[i] <<endl;
    }

    cout << "soluzione sistema B:" <<endl;
    for(int i = 0; i<nb; i++){
        cout << xb[i] <<endl;
    }

    cout << "soluzione sistema C:" <<endl;
    for(int i = 0; i<nb; i++){
        cout << xc[i] <<endl;
    }

    return 0;
}


void printMatrix(double **M, int nrow, int ncol) {
    cout << fixed << setprecision(4);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            cout << setw(10) << M[i][j] << "  ";
        }
        cout << endl;
    }
}


void eq_resolver(double **A, double * b, double *x, int n){

    for (int k = 0; k<n-1;k++) {            
        for (int i = k+1 ; i<n ; i++) {        
            double g = A[i][k]/A[k][k];
            for (int j = k+1 ; j<n ; j++) A[i][j] -= g*A[k][j];
            A[i][k] = 0.0;
            b[i]   -= g*b[k];
        }
    }  
    

    for (int i = n-1;i>= 0; i--){    
        double tmp = b[i];     
        for (int j = n-1;j > i; j--) tmp-= x[j]*A[i][j];    
        x[i] = tmp/A[i][i];  
    }

}



void eq_resolver_pp(double **A, double *b, double *x, int n) {

    // Fase di eliminazione
    for (int k = 0; k < n - 1; k++) {

        // --- PARTIAL PIVOTING ---
        int pivot_row = k;
        double max_val = fabs(A[k][k]);

        for (int i = k + 1; i < n; i++) {
            if (fabs(A[i][k]) > max_val) {
                max_val = fabs(A[i][k]);
                pivot_row = i;
            }
        }

        // Se necessario, scambiamo le righe
        if (pivot_row != k) {
            double *tmp_row = A[k];
            A[k] = A[pivot_row];
            A[pivot_row] = tmp_row;

            double tmp_b = b[k];
            b[k] = b[pivot_row];
            b[pivot_row] = tmp_b;
        }
        // -------------------------

        // Eliminazione in stile Gauss
        for (int i = k + 1; i < n; i++) {
            double g = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++)
                A[i][j] -= g * A[k][j];

            A[i][k] = 0.0;
            b[i] -= g * b[k];
        }
    }

    // Fase di sostituzione all'indietro
    for (int i = n - 1; i >= 0; i--) {
        double tmp = b[i];
        for (int j = n - 1; j > i; j--)
            tmp -= A[i][j] * x[j];

        x[i] = tmp / A[i][i];
    }
}
