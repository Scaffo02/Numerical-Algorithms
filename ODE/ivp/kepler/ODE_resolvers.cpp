
/*==================================================================================================

                                   ODE_RESOLVERS.CPP   

Descrizione:
---------------------------------------------------------------------------------
Integrazione numerica delle equazioni del moto per una particella in un campo
gravitazionale centrale newtoniano, in forma adimensionale.

Le quantità fisiche (in unità cgs) vengono scalate tramite:
    L_ref = 1.5e13 cm   (≈ 1 AU)
    GM    = 1.327e26 cm³/s²
    T_ref = sqrt(L_ref³ / GM)

I risultati (tempo, posizione e velocità) vengono riconvertiti in unità fisiche
e salvati nel file "Traiettoria_dim.dat".

Metodi disponibili:
    1) Eulero
    2) Runge-Kutta 2
    3) Runge-Kutta 4

---------------------------------------------------------------------------------
Autore:
    Scaffardi Alessandro
Data:
    Ottobre 2025

====================================================================================================*/




#include "header.h"


//------------------------------------- ODE_RESOLVERS ----------------------------------------
void ODE_resolvers(int neq, double *Y, double dt, double t, int Nsteps) {
    ofstream file("Traiettoria.dat");

    int metodo;
     cout <<"Scegli il metodo da usare: \n "<<
                                        " 1) Eulero \n "<<
                                        " 2) RK_2   \n "<< 
                                        " 3) RK_4   \n "<<
                                        " 4) Position Velet \n "<<
                                        " 5) Velocity Velet "<<endl;
    cin >> metodo;

    double x[neq];  
    double v[neq];

    for(int i = 0; i < (neq); i++){

            x[i] = Y[i];
            v[i] = Y[i+neq];
    }

    for(int i=0; i<Nsteps; i++) {
        // Scrittura delle quantità dimensionali
        file << t*T_ref << "\t"          // tempo dimensionale [s]
             << Y[0]*L_ref << "\t"       // x dimensionale [cm]
             << Y[1]*L_ref << "\t"       // y dimensionale [cm]
             << Y[2]*L_ref/T_ref << "\t" // vx dimensionale [cm/s]
             << Y[3]*L_ref/T_ref << endl;// vy dimensionale [cm/s]

        // Aggiornamento secondo il metodo scelto
        switch(metodo) {
            case 1: EuleroStep(t, Y, dYdt, dt, neq); break;
            case 2: RK2(t, Y, dYdt, dt, neq); break;
            case 3: RK4(t, Y, dYdt, dt, neq); break;
            case 4: 
                    position_velet(x,v,dt,neq); 
                    for(int i = 0; i < (neq); i++){
                        Y[i] = x[i];
                        Y[i+neq] = v[i];
                    }
                    break;
            case 5: 
                    velocity_velet(x,v,dt,neq); 
                    for(int i = 0; i < (neq); i++){
                        Y[i] = x[i];
                        Y[i+neq] = v[i];
                    }
                    break;

            default: cout << "Metodo non valido!" << endl; return;
        }

        //double theta = 0.08;
        //dt = theta* (Y[0]*Y[0]+Y[1]*Y[1])/(Y[2]*Y[2]+Y[3]*Y[3]);
        t +=dt;
        
    }

    file.close();
}







//------------------------------------ METODO DI EULERO ----------------------------------------------

void EuleroStep(double t, double *Y, void (*RSHFunc)(double, double*, double*), double dt, int neq){
/*
  Esegue un singolo passo di integrazione con il metodo di Eulero per il sistema dY/dt = RHS.
 
  - t         Tempo corrente.
  - Y         Vettore delle variabili del sistema (modificato in uscita).
  - RHSFunc   Puntatore alla funzione che calcola il termine destro del sistema (es. dYdt()).
  - dt        Passo temporale di integrazione.
  - neq       Numero di equazioni del sistema (dimensione del vettore Y[]).
*/
    int k = 0;
    double rhs[256];

    RSHFunc(t, Y, rhs);
    for(k=0; k<neq; k++){
        Y[k] += dt*rhs[k];
    }
}



//------------------------------------------ RK2 ------------------------------------------------

void RK2(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){
/**
  Esegue un singolo passo di integrazione con il metodo di Runge-Kutta del 2° ordine (RK2)
  per il sistema dY/dt = RHS.
 
  - t         Tempo corrente.
  - Y         Vettore delle variabili del sistema (modificato in uscita).
  - RHSFunc   Puntatore alla funzione che calcola il termine destro del sistema (es. dYdt()).
  - dt        Passo temporale di integrazione.
  - neq       Numero di equazioni del sistema (dimensione del vettore Y[]).
 */
    double Y1[neq], k1[neq], k2[neq];

    RHSFunc(t,Y,k1);
    for(int i = 0; i<neq ; i++){
        Y1[i] = Y[i] + 0.5*dt*k1[i]; 
    }

    RHSFunc(t+0.5*dt,Y1,k2);
    for(int i = 0; i<neq ; i++){
        Y[i] += dt*k2[i];
    }

}



//------------------------------------------ RK4 ------------------------------------------------
/**
  Esegue un singolo passo di integrazione con il metodo di Runge-Kutta del 4° ordine (RK2)
  per il sistema dY/dt = RHS.
 
  - t         Tempo corrente.
  - Y         Vettore delle variabili del sistema (modificato in uscita).
  - RHSFunc   Puntatore alla funzione che calcola il termine destro del sistema (es. dYdt()).
  - dt        Passo temporale di integrazione.
  - neq       Numero di equazioni del sistema (dimensione del vettore Y[]).
 */
void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){

    double Y1[neq], Y2[neq], Y3[neq], k1[neq], k2[neq], k3[neq], k4[neq];

    RHSFunc(t,Y,k1);
    for(int i = 0; i<neq ; i++){
        Y1[i] = Y[i] + dt*k1[i]/2.0;
    }

    RHSFunc(t+dt/2,Y1,k2);
    for(int i = 0; i<neq ; i++){
        Y2[i] = Y[i] + dt*k2[i]/2.0;
    }

    RHSFunc(t+dt/2,Y2,k3);
    for(int i = 0; i<neq ; i++){
        Y3[i] = Y[i] + dt*k3[i];
    }
    
    RHSFunc(t+dt,Y3,k4);
    for(int i = 0; i<neq ; i++){
        Y[i] = Y[i] + (dt/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    }


}

//------------------------------- POSITION VELET ------------------------------

void position_velet(double *x, double *v, double dt,int neq){

    double a[neq];

    

    for(int i = 0; i < neq ; i++){
        x[i] += 0.5*dt*v[i];
    }

    acc(x,a,neq);

    for(int i = 0; i < neq ; i++){
        v[i] += dt*a[i];
    }

    for(int i = 0; i < neq ; i++){
        x[i] += 0.5*dt*v[i];
    }
    
}



//----------------------------- VELOCITY VELET --------------------------------

void velocity_velet(double *x, double *v, double dt,int neq){

    double a[neq];

    acc(x,a,neq);

    for(int i = 0; i < neq ; i++){
        v[i] += 0.5*dt*a[i];
    }

    for(int i = 0; i < neq ; i++){
        x[i] += dt*v[i];
    }

    acc(x,a,neq);
    for(int i = 0; i < neq ; i++){
        v[i] += 0.5*dt*a[i];
    }

    
}



//---------------------------- ACCELERATION -----------------------------------

void acc(double *x,double *a,int neq){

    double r = sqrt(x[0]*x[0]+x[1]*x[1]);
    double factor = GM * T_ref*T_ref / (L_ref*L_ref*L_ref);

    for(int i = 0; i<(neq); i++){
        a[i] = -factor*x[i]/(r*r*r);
    }
    
    
}


//------------------------ FUNZIONE RHS (dimensionless) --------------------------------
void dYdt(double t, double *Y, double *k) {
    double x = Y[0];
    double y = Y[1];
    double vx = Y[2];
    double vy = Y[3];
    
    double r = sqrt(x*x + y*y);

    // fattore per passare da dimensionale a dimensionless
    double factor = GM * T_ref*T_ref / (L_ref*L_ref*L_ref);

    k[0] = vx;
    k[1] = vy;
    k[2] = - factor * x / (r*r*r);
    k[3] = - factor * y / (r*r*r);
}






