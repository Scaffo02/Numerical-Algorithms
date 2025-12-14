#include "header_integral.h"

int main(){

    double a = 0;
    double b = 0;
    int N = 0;
    int Ng = 3;
    int i = 0;

    cout <<"Inserisci estremo di integrazione sinistro:" <<endl;
    cin >> a;

    cout <<"Inserisci estremo di integrazione destro:" << endl;
    cin >> b;

    if(b<a){
        cout <<"L'estremo destro è più piccolo del sinistro!!!"<<endl;
        return 0;
    }

    N = abs( int(b-a)*10 );

    cout <<"Quale metodo di integrazione vuoi usare??\n 1) Rettangoli \n 2) Trapezi \n 3) Simpson \n 4) Gauss"<<endl;
    cin >> i;

    switch(i){

    case 1:
        cout <<"Risultato con metodo dei rettangoli: "<< rettangoli(f,a,b,N)<<endl;
        break;
    case 2:
        cout <<"Risultato con metodo dei trapezi: "<< trapezi(f,a,b,N)<<endl;
        break;
    case 3:
        cout <<"Risultato con metodo dei simpson: "<< simpson(f,a,b,N)<<endl;
        break;
    case 4:
        cout <<"Risultato con metodo dei gauss: "<< gauss(f,a,b,N,Ng)<<endl;
        break;

    default:
        cout <<"Numero non valido!!!"<<endl;
        return 0;

    }

    return 0;

}
