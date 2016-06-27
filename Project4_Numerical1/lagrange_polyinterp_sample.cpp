#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
      
using namespace std;
      
      double fx_smooth(double x) {
      return 1.0/(1.0+x*x);
      }


       // return y=P_n(x)   P_n interpolates the points in xsten.
      void lagrange_interp(double* xsten,double (*fx)(double),
       int n,double x,double& y) {

      double* fsten=new double[n+1];
      int i,j;
      double L;

      for (i=0;i<=n;i++) {
       fsten[i]=fx(xsten[i]);
      }
      y=0.0;
      for (i=0;i<=n;i++) {
       L=1.0;
       for (j=0;j<=n;j++) {
        if (i!=j) {
         L=L*(x-xsten[j])/(xsten[i]-xsten[j]);
        }
       }
       y=y+fsten[i]*L;
      }

      delete fsten;

      }
   
      int main() {

       int n=10;

       double xlo=-5.0;
       double xhi=5.0;
       double h=(xhi-xlo)/n;
       double* xsten=new double[n+1];
       for (int i=0;i<=n;i++) {
        xsten[i]=xlo+i*h;
       }

       ofstream lagrangefile;
       lagrangefile.open("lagrange.dat");
       ofstream originalfile;
       originalfile.open("original.dat");
       cout << "output filenames are lagrange.dat and original.dat \n";
  
       int nplot=1000;
       double hplot=(xhi-xlo)/nplot;

       for (int i=0;i<=nplot;i++) {
        double xplot=xlo+i*hplot;  //xplot = xlo + i*(xhi-xlo)/nplot
        double fplot;
        double foriginal;
        lagrange_interp(xsten,fx_smooth,n,xplot,fplot);
        lagrangefile << xplot << ' ' << fplot << '\n';
        foriginal=fx_smooth(xplot);
        originalfile << xplot << ' ' << foriginal << '\n';
       }

       lagrangefile.close(); 
       originalfile.close(); 
  
       delete xsten;

      } 
