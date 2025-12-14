#include <iostream>
#include <iomanip>

using namespace std;
void matmul(double **A, double **B, double **C,
            int n, int m, int p);
void printMatrix(double **M, int nrow, int ncol);

int main() {
    
    int n = 4, m = 4, p = 4;

    // Alloca matrici A, B e C
    double **A = new double*[n];
    double **B = new double*[m];
    double **C = new double*[n];

    A[0] = new double[n*m];
    B[0] = new double[m*p];
    C[0] = new double[n*p];

    for(int i = 1; i < n; i++) A[i] = A[i-1] + m;
    for(int i = 1; i < m; i++) B[i] = B[i-1] + p;
    for(int i = 1; i < n; i++) C[i] = C[i-1] + p;

    // inizializza A e B
    for(int i = 0; i < n; i++)
        for(int j = 0; j < m; j++)
            A[i][j] = i + j;

    for(int i = 0; i < m; i++)
        for(int j = 0; j < p; j++)
            B[i][j] = (i == j ? 1 : 0); // matrice identità

    // prodotto
    matmul(A, B, C, n, m, p);

    cout << "Risultato C = A * B:" << endl;
    printMatrix(C, n, p);

    return 0;
}


// Moltiplica C = A * B
// A: n x m
// B: m x p
// C: n x p  (deve essere già allocata)
void matmul(double **A, double **B, double **C,
            int n, int m, int p) {

    // inizializza C a zero
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            C[i][j] = 0.0;

    // prodotto matriciale
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < m; k++) {
            double aik = A[i][k];
            for (int j = 0; j < p; j++) {
                C[i][j] += aik * B[k][j];
            }
        }
    }
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
