/*
Program that solves the two-dimensional Ising model using the Metropolis algorithm.
There is no external magnetic field and the coupling constant J is set to zero.
The Boltzmann constant is set to zero, so temperature has dimension energy.*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

// output file
ofstream ofile;

//Function to initialize spin matrix and compute initial energy and magnetization
void initialize(int, int **, double&, double&);
// Function that ensures periodic boundary conditions
int periodic(int, int);
// Function for the Metropolis algo
void metropolis(int, int, double, double *);
// prints results to file
void WriteToFile(int, int, double, double *);


int main(int argc, char* argv[]){
    string filename;
    int n_spins, MCcycles;
    double T_min, T_max, T_step;
    filename = argv[1];
    n_spins = atoi(argv[2]);
    MCcycles = atoi(argv[3]);
    double expectation_values[5];    // array to hold expectation values of E, E^2, M, M^2 and |M|
    double T = 1.0;

    metropolis(n_spins, MCcycles, T, expectation_values);
    WriteToFile(n_spins, MCcycles, T, expectation_values);

}

//Function that sets up spin matrix and computes initial energy and magnetization
void initialize(int n_spins, int ** spins, double& E, double& M){

    // Initializing spin matrix with all spin up and computing magnetization
    for (int i = 0; i < n_spins; i++){
        for (int j = 0; j < n_spins; j++){
            spins[i][j] = 1.0;
            M += (double) spins[i][j];
        }
    }
    // compute initial energy
    for (int i = 0; i < n_spins; i++){
        for (int j = 0; j < n_spins; j++){
            E -= (double) (spins[i][j]*(spins[i][periodic(j+1, n_spins)]+spins[i][j]*spins[periodic(i+1,n_spins)][j]));
        }
    }
}

void metropolis(int n_spins, int MCcycles, double temp, double *expectation_values){
    int ** spins;                    // spin matrix
    int deltaE[5];                   // array to store possible energy differences
    double exp_deltaE[5];            // array to store precalculated w = exp(-dE/T)
    double T = temp;                 // temperature, kT/J

    random_device rd;     // random device
    mt19937_64 rng(rd());
    uniform_real_distribution<double> uniform(0.0, 1.0);  //Uniform distribution of number between 0 and 1

    // Allocate memory for 2D spin matrix
    spins = new int * [n_spins];
    for (int i = 0; i < n_spins; i++){
       spins[i] = new int [n_spins];
    }
    // Set energy and magnetization to zero
    double E = 0.0;
    double M = 0.0;

    // Set up spin matrix and compute initial energy and magnetization
    initialize(n_spins, spins, E, M);

    // Define possible energy differences
    deltaE[0] = -16; deltaE[1] = -8; deltaE[2] = 0; deltaE[3] = 8; deltaE[4] = 16;

    // Set up array with the exponential of all the different energy differences
    for (int i=0; i<5; i++){
        exp_deltaE[i] = exp(-(double)deltaE[i]/T);
    }

    // start Monte Carlo cycles
    for (int cycles = 1; cycles <= MCcycles; cycles ++){

        // Metropolis algorithm
        for (int i = 0; i < n_spins; i++){
            for (int j = 0; j < n_spins; j++){
                // Choose random matrix element (i_index,j_index)
                int i_index = (int) (uniform(rng)*(double)n_spins);
                int j_index = (int) (uniform(rng)*(double)n_spins);
                // Compute difference in energy after flipping spin at (i_index,j_index)
                int dE = 2*spins[i_index][j_index]*(spins[i_index][periodic(j_index+1, n_spins)]+spins[periodic(i_index+1, n_spins)][j_index]
                        +spins[i_index][periodic(i_index-1, n_spins)]+spins[periodic(i_index-1, n_spins)][j_index]);


                if (dE <= 0){
                    spins[i_index][j_index] *= -1;
                    E += (double) dE;
                    M += (double) (2*spins[i_index][j_index]);
                }
                else{
                    double w = 0;
                    for (int i = 0; i < 5; i++){
                        if (dE == deltaE[i]){
                           w = exp_deltaE[i];             // Probability for transition between states
                        }
                    }
                    double rand = (double)(uniform(rng));
                    if (rand <= w){
                        spins[i_index][j_index] *= -1;
                        E += (double) dE;
                        M += (double) (2*spins[i_index][j_index]);
                    }
                }
            }
        }
        // Update expectation values
        expectation_values[0] += E;
        expectation_values[1] += E*E;
        expectation_values[2] += M;
        expectation_values[3] += M*M;
        expectation_values[4] += fabs(M);


    }
}

// Function that ensures periodic boundary conditions by moving to the first element of the matrix if the index is out of bounds
int periodic(int index, int n_spins){
    if (index >= n_spins){
        return 0;
    }
    else if (index == -1){
        return n_spins - 1;
    }
    else {return index;}
}

// Function that writes results to file
void WriteToFile(int n_spins, int MCcycles, double temp, double * expectation_values){
    double norm = 1.0/(double)(MCcycles);      // divide by number of cycles
    double E_expected = expectation_values[0]*norm;
    double E2_expected = expectation_values[1]*norm;
    double M_expected = expectation_values[2]*norm;
    double M2_expected = expectation_values[3]*norm;
    double Mabs_expected = expectation_values[4]*norm;

}




