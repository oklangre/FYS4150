#include <iostream>
#include <cstdlib>

double compute_energy(int dim, double ** A, double J);
using namespace std;

int main(){
    int L = 2;       // number of spins
    double T = 1.0;  // kT/J
    double J = 1.0;


    // setting up emtpy spin matrix, to be improved
    double ** A = new double * [L];
    for (int i = 0; i < L; i++){
       A[i] = new double [L];
    }
    // Initializing spin matrix with all spin up
    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            A[i][j] = 1.0;
        }
    }

    int sweeps = 0;
    for (int i=0; i<50; i++){
       int num = rand()%2;
       int num2 = rand()%2;
       A[num][num2] = -1.0*A[num][num2];
       compute_energy(L, A, J);



    }




}

double compute_energy(int dim, double ** A, double J){
    int im = dim-1;
    int jm = dim-1;
    double E = 0;
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            E -= J*A[i][j]*A[i][jm] + J*A[i][j]*A[im][j];
            jm = j;
            cout << E << endl;
        }
        im = i;
    }
    return E;
}

