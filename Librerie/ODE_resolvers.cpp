
/*==================================================================================================

                                         ODE_RESOLVERS.CPP

Questa libreria implementa diversi metodi numerici per la risoluzione di sistemi di equazioni
differenziali ordinarie (ODE) a partire da condizioni iniziali assegnate.

  Sono disponibili vari schemi di integrazione, tra cui:
  - Metodo di Eulero
  - Metodo di Runge-Kutta del 2° ordine (RK2)
  - Metodo di Runge-Kutta del 4° ordine (RK4)
  - Position Verlet
  - Velocity Verlet

La funzione principale `ODE_resolvers` seleziona e applica automaticamente il metodo desiderato
per risolvere il sistema di ODE specificato.

Per utilizzare correttamente questa libreria, è necessario:
1. Definire nella funzione `dY/dt` il sistema di equazioni differenziali da risolvere.
2. Definire nella funzione 'acc' il termine di accelerazione
2. Includere il file nel proprio header tramite:
   `#include "MY_HEADER"`

**Attenzione**
Assicurarsi che la funzione `dY/dt` sia coerente con la dimensione del sistema e che le condizioni
iniziali siano definite correttamente. Idem per la funzione 'acc'.

Autore:
Alessandro Scaffardi

==================================================================================================*/


#include "header.h"



//------------------------------------- ODE_RESOLVERS ----------------------------------------

void ODE_resolvers(int neq, double t, double dt, double *Y,int Nsteps){
    ofstream file("Soluzione.dat");

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
        // Scrittura delle quantità su file
        file << t << "\t"          // tempo 
             << Y[0] << "\t"       // x 
             << Y[1] << "\t"       // y 
             << Y[2] << "\t"       // vx
             << Y[3] << endl;      // vy

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

void RK4(double t, double *Y, void (*RHSFunc)(double, double*, double*), double dt, int neq){
/**
  Esegue un singolo passo di integrazione con il metodo di Runge-Kutta del 4° ordine (RK2)
  per il sistema dY/dt = RHS.
 
  - t         Tempo corrente.
  - Y         Vettore delle variabili del sistema (modificato in uscita).
  - RHSFunc   Puntatore alla funzione che calcola il termine destro del sistema (es. dYdt()).
  - dt        Passo temporale di integrazione.
  - neq       Numero di equazioni del sistema (dimensione del vettore Y[]).
 */
    double Y1[neq], Y2[neq], Y3[neq], k1[neq], k2[neq], k3[neq], k4[neq];

    RHSFunc(t,Y,k1);
    for(int i = 0; i<neq ; i++){
        Y1[i] = Y[i] + dt*k1[i]*0.5;
    }

    RHSFunc(t+dt/2,Y1,k2);
    for(int i = 0; i<neq ; i++){
        Y2[i] = Y[i] + dt*k2[i]*0.5;
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
/**
  Aggiorna le posizioni e le velocità di un sistema dinamico utilizzando lo 
  schema di integrazione di Position Verlet.
  - x    Vettore delle posizioni del sistema.
  - v    Vettore delle velocità correnti del sistema.
  - dt   Passo temporale di integrazione.
  - neq  Numero di equazioni/particelle (dimensione dei vettori x[] e v[]).
 */
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
/**
  Aggiorna le posizioni e le velocità di un sistema dinamico utilizzando lo 
  schema di integrazione di Velocity Verlet.

  - x    Vettore delle posizioni correnti del sistema.
  - v    Vettore delle velocità del sistema.
  - dt   Passo temporale di integrazione.
  - neq  Numero di equazioni/particelle (dimensione dei vettori x[] e v[]).
 */
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
/**
  Calcola le accelerazioni del sistema a partire dalle posizioni correnti,
  tipicamente applicando le leggi dinamiche o le forze definite dal modello fisico.

  - x    Vettore delle posizioni correnti del sistema.
  - a    Vettore delle accelerazioni.
  - neq  Numero di equazioni/particelle (dimensione dei vettori x[] e a[]).
 */
void acc(double *x,double *a,int neq){

    for(int i = 0; i<neq; i++){
        
        a[i] = -omega*omega*x[i];

    }
    
}


//-------------------------------------- ODE ----------------------------------
/*
 Calcola il termine destro del sistema di equazioni differenziali dY/dt = k.
 
 - t   Tempo corrente.
 - Y   Vettore delle variabili del sistema (Y[0] = x, Y[1] = y).
 - k   Vettore risultante contenente le derivate (k[0] = dx/dt, k[1] = dy/dt).
 
  In questo esempio:
     dx/dt = y
     dy/dt = -x
 */
void dYdt(double t, double *Y, double *k){

    double x = Y[0], y = Y[1];
    k[0] = y;
    k[1] = -x;
    
}



