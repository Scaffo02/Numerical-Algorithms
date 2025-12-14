#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector> // Necessario per salvare le soluzioni se volessimo, o per logica vettoriale

using namespace std;

// --- Parametri globali ---
const double v0 = 200.0;     // velocità iniziale [m/s]
const double target_x = 1000.0; // distanza bersaglio [m]
const double g = 9.81;
const double B = 0.0005;     // Coefficiente attrito (normalizzato dalla massa)
const double Delta_t = 0.0005;
const int NEQ = 4;           // Numero di equazioni ODE

// Prototipi funzioni
void dYdt(double t, const double *Y, double *k);
void RK4(double t, double *Y, double dt);
double Range(double angle_rad); 
double Residual(double angle_rad);
double Bisection(double (*func)(double), double a, double b, double tol);
void PrintTrajectory(double angle_rad, const char* filename);
void BracketingAndSolve(double start, double end, int N_subdivisions, double tol);
double FindMaxRangeAngle(double a_low, double a_high, double Delta_angle);
void GenerateRangePlotData(double a_low, double a_high, double Delta_angle); 


// ====================================== MAIN =================================
int main(){
    double tol = 1e-6;
    
    // --- 1. Ricerca Automatica delle Soluzioni (Bracketing + Bisezione) ---
    cout << "--- Avvio Bracketing e Ricerca Radici ---" << endl;
    cout << "Intervallo di ricerca: 0 a 90 gradi" << endl;
    cout << "Bersaglio a: " << target_x << " m\n" << endl;

    // Cerchiamo radici tra 0 e PI/2 (0 e 90 gradi) dividendo in 100 intervalli
    // Se ci sono soluzioni, la funzione le troverà, le stamperà e salverà i file.
    BracketingAndSolve(0.0, M_PI / 2.0, 100, tol);

    // --- 2. Ricerca Angolo di Gittata Massima (Opzionale, utile per analisi) ---
    cout << "\n--- Analisi Gittata Massima ---" << endl;
    double theta_max = FindMaxRangeAngle(0.0, M_PI / 2.0, 0.0001); 

    if (theta_max > 0.0) {
        double max_range = Range(theta_max);
        
        cout << fixed << setprecision(4);
        cout << "Angolo per gittata massima: " << theta_max * 180.0 / M_PI << " gradi" << endl;
        cout << "Gittata Massima raggiungibile: " << max_range << " metri" << endl;
        
        if (max_range < target_x) {
            cout << "ATTENZIONE: Il bersaglio e' fuori dalla portata massima!" << endl;
        } else {
            PrintTrajectory(theta_max, "traj_max_gittata.dat");
        }
    }

    // --- 3. Generazione Dati per Grafico Gittata vs Angolo ---
    GenerateRangePlotData(0.0, M_PI / 2.0, 0.005); 
    cout << "\nDati per il grafico della gittata salvati in gittata_vs_angolo.dat" << endl;

    // --- 4. Salvataggio Bersaglio ---
    ofstream ftarget("target.dat");
    if (ftarget.is_open()) {
        ftarget << target_x << " " << 0.0 << endl;
        ftarget.close();
    }

    return 0;
}



// =================================== IMPLEMENTAZIONE FUNZIONI ================

// ---------------- FUNZIONE DI BRACKETING ------------------------------
/*
  Questa funzione divide l'intervallo totale [start, end] in N sotto-intervalli.
  Controlla se c'è un cambio di segno del Residuo agli estremi di ogni sotto-intervallo.
  Se c'è cambio di segno -> Chiama la Bisezione.
*/
void BracketingAndSolve(double start, double end, int N_subdivisions, double tol) {
    double dx = (end - start) / N_subdivisions;
    double x_L, x_R;
    int solution_count = 0;

    for (int i = 0; i < N_subdivisions; i++) {
        x_L = start + i * dx;
        x_R = start + (i + 1) * dx;

        double y_L = Residual(x_L);
        double y_R = Residual(x_R);

        // Controllo cambio di segno (Teorema degli Zeri)
        if (y_L * y_R <= 0.0) {
            solution_count++;
            
            // Determina l'etichetta (solo estetica)
            string label;
            string filename;
            
            // Logica semplice: se l'angolo è basso (< 45 gradi circa) è tiro teso
            // Nota: con l'attrito l'angolo di max gittata è < 45, ma usiamo 45 come riferimento grossolano
            double mid_deg = ((x_L + x_R) / 2.0) * (180.0 / M_PI);
            
            if (mid_deg < 45.0) {
                label = "Soluzione Bassa (Teso)";
                filename = "traj_bassa.dat";
            } else {
                label = "Soluzione Alta (Mortaio)";
                filename = "traj_alta.dat";
            }

            cout << "-> Radice trovata nell'intervallo [" 
                 << x_L * 180/M_PI << ", " << x_R * 180/M_PI << "] gradi." << endl;

            // --- CHIAMATA ALLA BISEZIONE ---
            double root = Bisection(Residual, x_L, x_R, tol);

            cout << "   " << label << ": " << root * 180.0 / M_PI << " gradi." << endl;
            cout << "   Errore residuo: " << Residual(root) << " m" << endl;
            
            PrintTrajectory(root, filename.c_str());
            cout << "--------------------------------------------" << endl;
        }
    }

    if (solution_count == 0) {
        cout << "Nessuna soluzione trovata nell'intervallo specificato." << endl;
        cout << "Il bersaglio potrebbe essere troppo lontano." << endl;
    }
}

//-------------------------------------- ODE ----------------------------------
void dYdt(double t, const double *Y, double *k){
    double vx = Y[2];
    double vy = Y[3];
    double v = sqrt(vx*vx + vy*vy);

    k[0] = vx;
    k[1] = vy;
    k[2] = -B * v * vx; 
    k[3] = -g - B * v * vy;
}

//------------------------------------------ RK4 -------------------------------
void RK4(double t, double *Y, double dt){
    double Y1[NEQ], Y2[NEQ], Y3[NEQ];
    double k1[NEQ], k2[NEQ], k3[NEQ], k4[NEQ];

    dYdt(t, Y, k1);
    for(int i = 0; i < NEQ; i++) Y1[i] = Y[i] + dt*k1[i]*0.5;

    dYdt(t+dt/2.0, Y1, k2);
    for(int i = 0; i < NEQ; i++) Y2[i] = Y[i] + dt*k2[i]*0.5;

    dYdt(t+dt/2.0, Y2, k3);
    for(int i = 0; i < NEQ; i++) Y3[i] = Y[i] + dt*k3[i];

    dYdt(t+dt, Y3, k4);
    for(int i = 0; i < NEQ; i++) {
        Y[i] = Y[i] + (dt/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    }
}

// ------------------------------------ GITTATA (Range) ------------------------
double Range(double angle_rad){
    double dt = Delta_t; 
    double t = 0;
    double Y[NEQ] = {0, 0, v0*cos(angle_rad), v0*sin(angle_rad)};

    while(Y[1] >= 0){
        RK4(t, Y, dt);
        t += dt;
        if (t > 150) return Y[0]; 
    }
    return Y[0];
}

//------------------------------------ RESIDUAL --------------------------------
double Residual(double angle_rad){
    return Range(angle_rad) - target_x;
}

//-------------------------- BISEZIONE (Root Finding) --------------------------
double Bisection(double (*func)(double), double a, double b, double tol){
    double fa = func(a);
    double xm = (a+b)/2;
    double fm = func(xm);
    int iter = 0;
    
    while(fabs(b-a) > tol && iter < 1000){
        xm = (a+b)/2;
        fm = func(xm);
        
        if(fm == 0.0) break;

        if(fa * fm < 0){
            b = xm;
        } else {
            a = xm;
            fa = fm; 
        }
        iter++;
    }
    return xm;
}

// -------------------------- RICERCA ANGOLO DI GITTATA MASSIMA ----------------
double FindMaxRangeAngle(double a_low, double a_high, double Delta_angle){
    double max_range = 0.0;
    double theta_max = -1.0;
    for (double angle = a_low; angle <= a_high; angle += Delta_angle) {
        double current_range = Range(angle);
        if (current_range > max_range) {
            max_range = current_range;
            theta_max = angle;
        }
    }
    return theta_max;
}

// -------------------------- GENERAZIONE DATI PER IL PLOT ---------------------
void GenerateRangePlotData(double a_low, double a_high, double Delta_angle){
    ofstream fout("gittata_vs_angolo.dat");
    if (!fout.is_open()) return;
    
    fout << fixed << setprecision(4);
    fout << "# Angolo [Gradi] Gittata [m]" << endl;
    
    for (double angle_rad = a_low; angle_rad <= a_high; angle_rad += Delta_angle) {
        fout << angle_rad * 180.0 / M_PI << " " << Range(angle_rad) << "\n";
    }
    fout.close();
}

//-------------------------------------- TRAIETTORIA ---------------------------
void PrintTrajectory(double angle_rad, const char* filename){
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Errore file: " << filename << endl;
        return;
    }
    
    double dt = Delta_t; 
    double t = 0;
    double Y[NEQ] = {0.0, 0.0, v0*cos(angle_rad), v0*sin(angle_rad)};

    fout << fixed << setprecision(4) << "# x[m] y[m]" << endl;
    while(Y[1] >= 0){
        fout << Y[0] << " " << Y[1] << "\n";
        RK4(t, Y, dt);
        t += dt;
    }
    cout << "Traiettoria salvata in " << filename << endl;
    fout.close();
}