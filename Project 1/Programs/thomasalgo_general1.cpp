// Program to solve a matrix equation of form Av=y for any value of the matrix elements

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;
ofstream ofile;

// Source term  
inline double f(double x){return 100.0*exp(-10.0*x);}

// Exact solution for v
inline double u(double x){return 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x);}

int main(int argc, char* argv[]){
   if (argc <= 5){
   cout << "Bad usage: write also filename for output file, number of gridpoints, value of element on left diagonal, value of element on middle diagonal and value of element on right diagonal." << endl;
   exit(1);
   }  
   int n = atoi(argv[2]);
   string filename = argv[1];
   // Value of matrix elements from terminal
   double a_value = atof(argv[3]);
   double b_value = atof(argv[4]);    
   double c_value = atof(argv[5]);
   
   double* b_tilde = new double [n+1]; double* y = new double[n+1]; double* y_tilde = new double[n+1];
   double* x = new double[n+1]; double* v = new double[n+1];
   double* a = new double[n]; double* b = new double[n+1]; double* c = new double[n];
   string fileout = filename;

   // Defining diagonals
   for (int i = 0; i < n-1; i++){
      a[i] = a_value;
      c[i] = c_value;
   }
   for (int i=0; i<n; i++){
      b[i] = b_value;
   }
   
   
   double h = 1.0/n;
   double hh = h*h;
 
   // Calculating d(x) = h^2f(x)
   for (int i=0; i<=n; i++){
      x[i] = i*h; 
      y[i] = hh*f(i*h);}
   // Imposing boundary conditions
   v[0] = 0.0; v[n] = 0.0;
   // Defining first and last elements of updated vectors
   b_tilde[0] = b[0];
   b_tilde[n] = b[n];
   y_tilde[0] = y[0];
   y_tilde[n] = y[n];

   // Forward substitution and updating of coefficients 
   for (int i = 1; i < n; i++){
      b_tilde[i] = b[i] - a[i-1]*c[i-1]/b_tilde[i-1];
      y_tilde[i] = y[i] - a[i-1]*y_tilde[i-1]/b_tilde[i-1];
   }
   // Backward substitution to solve for v_i
   v[n-1] = y_tilde[n-1]/b_tilde[n-1];
   for (int i = n-2; i>0; i--){ v[i] = y_tilde[i]-c[i]*v[i+1]/b_tilde[i];}  

   ofile.open(fileout.c_str());
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
	double xval = x[i];
 	 double RelativeError = fabs((u(xval)-v[i])/u(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << v[i];
         ofile << setw(15) << setprecision(8) << u(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] x; delete [] y_tilde; delete [] a; delete [] b; delete [] c; delete [] b_tilde; delete [] v;
   
   return 0;
}
  
