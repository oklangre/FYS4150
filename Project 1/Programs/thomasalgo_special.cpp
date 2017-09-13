// Program to solve a matrix equation of form Av=y when the diagonal elements are -1, 2 and 1.

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
#include <time.h>

using namespace std;
ofstream ofile;

// Source term  
inline double f(double x){return 100.0*exp(-10.0*x);}

// Exact solution for v
inline double u(double x){return 1.0-(1.0-exp(-10.0))*x-exp(-10.0*x);}

int main(int argc, char* argv[]){
   if (argc <= 2){
   cout << "Bad usage: write also filename for output file and number of gridpoints." << endl;
   exit(1);
   }  
   
   int n = atoi(argv[2]);
  
   string filename = argv[1];
   double b_value = 2.0;    // Value of diagonal matrix elements
   double* b_tilde = new double [n+1]; double* y = new double[n+1]; double* y_tilde = new double[n+1];
   double* x = new double[n+1]; double* v = new double[n+1];
   
   string fileout = filename;
   
   
   double h = 1.0/n;
   double hh = h*h;
 
   // Calculating d(x) = h^2f(x)
   for (int i=0; i<=n; i++){
      x[i] = i*h; 
      y[i] = hh*f(i*h);}
   // Imposing boundary conditions
   v[0] = 0.0; v[n] = 0.0;
   // Defining first and last elements of updated vectors
   b_tilde[0] = b_tilde[n] = b_value;
   y_tilde[0] =  y[0];
   y_tilde[n] = y[n];

   // Forward substitution 
   for (int i = 1; i < n; i++){
      b_tilde[i] = (i+1.0)/((double)i);
      y_tilde[i] = y[i]+(i-1)*y_tilde[i-1]/(double(i));
   }
   // Backward substitution
   v[n-1] = y_tilde[n-1]/b_tilde[n-1];
   for (int i = n-2; i>0; i--){ v[i] = ((double)i)/(i+1)*(y_tilde[i]+v[i+1]);}  

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
      delete [] x; delete [] y_tilde; delete [] b_tilde; delete [] v;
      
   
   return 0;
}
  
   

 
  
  
 
