#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
void printMatrix(double **M, int nrow, int ncol);
bool solve_tridiagonal_full(double **A, double *d, double *x, int n);

int main(){

    int n = 5;
    double **A = new double*[n];
    double b[n]={1,0,3,1,0};
    double x[n];
    A[0] = new double[n*n];
 
    for(int i = 1; i < n; i++) A[i] = A[i-1] + n;
    
    A[0][0] = 2,  A[0][1] = 1, A[0][2] = 0,  A[0][3] = 0, A[0][4] = 0; 
    A[1][0] = 1,  A[1][1] = 2, A[1][2] = 1,  A[1][3] = 0, A[1][4] = 0;
    A[2][0] = 0,  A[2][1] = 1, A[2][2] = 2,  A[2][3] = 1, A[2][4] = 0;
    A[3][0] = 0,  A[3][1] = 0, A[3][2] = 1,  A[3][3] = 2, A[3][4] = 1;
    A[4][0] = 0,  A[4][1] = 0, A[4][2] = 0,  A[4][3] = 1, A[4][4] = 2;

    cout << "matrice A:" <<endl;
    printMatrix(A,n,n);

    solve_tridiagonal_full(A,b,x,n);
    cout << "soluzione sistema A:" <<endl;
    for(int i = 0; i<n; i++){
        cout << x[i] <<endl;
    }

    return 0;
}


bool solve_tridiagonal_full(double **A, double *d, double *x, int n)
{
    // Controlla che la matrice sia tridiagonale
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (abs(i - j) > 1 && A[i][j] != 0.0) {
                cerr << "Errore: la matrice non è tridiagonale.\n";
                return false;
            }
        }
    }

    // Estrai i vettori a, b, c
    double *a = new double[n];   // sotto-diagonale
    double *b = new double[n];   // diagonale
    double *c = new double[n];   // sopra-diagonale

    for (int i = 0; i < n; i++) {
        b[i] = A[i][i];
        a[i] = (i > 0) ? A[i][i - 1] : 0.0;
        c[i] = (i < n - 1) ? A[i][i + 1] : 0.0;
    }

    double *cp = new double[n];
    double *dp = new double[n];

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double denom = b[i] - a[i] * cp[i - 1];
        if (abs(denom) < 1e-14) {
            cerr << "Pivot nullo! Sistema instabile.\n";
            delete[] a; delete[] b; delete[] c;
            delete[] cp; delete[] dp;
            return false;
        }

        cp[i] = (i < n - 1) ? c[i] / denom : 0.0;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
    }

    x[n - 1] = dp[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    delete[] a; delete[] b; delete[] c;
    delete[] cp; delete[] dp;

    return true;
}



void printMatrix(double **M, int nrow, int ncol){
    cout << fixed << setprecision(4);
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            cout << setw(10) << M[i][j] << "  ";
        }
        cout << endl;
    }
}