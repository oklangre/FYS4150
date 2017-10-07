#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <fstream>

using namespace std;
ofstream ofile;
 
int main(int argc, char* argv[]){
   double x0 = 9.7760E-1;    // AU
   double y0 = 2.259788E-1;  // AU
   double vx0 = -4.0624E-3*365;  // AU/year
   double vy0 = 1.671862E-2*365; // AU/year
   double x1 = 9.7340257E-1;   // AU
   double y1 = 2.419834645E-1;

   double tmax = atof(argv[2]);   // years

   double t_step = 1/365.0;   // 1 day
   int n = tmax/t_step;
   cout << n << endl;
   string filename = argv[1];
   string fileout = filename;
   double * x = new double[n+1];
   double * y = new double[n+1];
   y[0] = y0;
   x[0] = x0;
   x[1] = x1;
   y[1] = y1;
   double * vx = new double[n+1];
   double * vy = new double[n+1];
   vx[0] = vx0;
   vy[0] = vy0;
   
   
   //double Msun = 2.0E30;  // sun mass
   //double G = 6.67E-11;   // gravitational constant, Nm^2/kg^2

   double pi = acos(-1.0);
   double fourpi2 = 4*pi*pi;

   ofile.open(fileout.c_str());
   ofile << setiosflags(ios::showpoint | ios:: uppercase);

   // Forward Euler
   for (int i = 0; i<n; i++){
      double r = sqrt(y[i]*y[i]+x[i]*x[i]);
      double ax = -fourpi2*x[i]/(r*r*r); 
      double ay = -fourpi2*y[i]/(r*r*r);
      vx[i+1] = vx[i]+ax*t_step;
      vy[i+1] = vy[i]+ay*t_step;
      x[i+1] = x[i]+vx[i+1]*t_step;
      y[i+1] = y[i] + vy[i+1]*t_step;
      ofile << setw(12) << x[i] << setw(20) << y[i] << endl;
   
   }
   ofile.close();
   
   /*
   // Velocity Verlet
   for (int i = 1; i<n; i++){
      double r = sqrt(y[i]*y[i]+x[i]*x[i]);
      double ax = -fourpi2*x[i]/(r*r*r); 
      double ay = -fourpi2*y[i]/(r*r*r);
      x[i+1] = 2*x[i]-x[i-1]+(t_step*t_step*ax);
      y[i+1] = 2*y[i]-y[i-1]+(t_step*t_step*ay);
      ofile << setw(12) << x[i] << setw(20) << y[i] << endl;
   
   }
   ofile.close();
   */

   
}   
   
   
    
   
   
