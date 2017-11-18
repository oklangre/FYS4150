#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

// output file
ofstream ofile;

void initialize(int, int **, double&, double&);
int periodic(int, int);

int main(int argc, char* argv[]){
    char *outfilename;
    int ** spins, n_spins;
    double E, M;
    n_spins = atoi(argv[1]);
    double expectation_values[5];    // array to hold expectation values of E, E^2, M, M^2 and |M|

    random_device rd;     // random device
    mt19937_64 rng(rd());
    uniform_real_distribution<double> uniform(0.0, 1.0);  //Uniform distribution of number between 0 and 1


    // Allocate memory for 2D spin matrix
    spins = new int * [n_spins];
    for (int i = 0; i < n_spins; i++){
       spins[i] = new int [n_spins];
    }

}

//Function that sets up spin matrix and computes initial energy and magnetization
void initialize(int n_spins, int ** spins, double& E, double& M){

    // Initializing spin matrix with all spin up and computing magnetization
    for (int i = 0; i < n_spins; i++){
        for (int j = 0; j < n_spins; j++){
            spins[i][j] = 1.0;
            M += (double) spins[i][j];
            cout << M << endl;
        }
    }
    // compute initial energy
    for (int i = 0; i < n_spins; i++){
        for (int j = 0; j < n_spins; j++){
            E -= (double) (spins[i][j]*(spins[i][periodic(j+1, n_spins)]+spins[i][j]*spins[periodic(i+1,n_spins)][j]));
        }
    }
}

// Function that ensures periodic boundary conditions by moving to the first element of the matrix if the index is out of bounds
int periodic(int index, int n_spins){
    if (index >= n_spins){
        return 0;
    }
    else {return index;}
}

void compute_energy_and_moment(int ** spins, int n_spins, double& E, double& M){
    for (int i = 0; i < n_spins; i++){
        for (int j = 0; j < n_spins; j++){
            M += (double) spins[i][j];
            E -= (double) (spins[i][j]*(spins[i][periodic(j+1, n_spins)]+spins[i][j]*spins[periodic(i+1,n_spins)][j]));
        }
    }

}
