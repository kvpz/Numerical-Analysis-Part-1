#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

      double fx_sample(double x) {
      return x*x;
      }
 
      void composite_midpoint(int n,double (*fx)(double), double a,double b,double& sum) {
      double xi,h;
      int i;

      std::cout << "calculating integral with a,b,n = " << a << ' ' <<
        b << ' ' << n << '\n';
      h=(b-a)/n;
      sum=0.0;
      for (i=0;i<n;i++) {
       xi=a+h*(i+0.5);
       sum=sum+fx(xi);
      }
      sum=h*sum;
     }
  
 
      int main() {

      int n;
      double a,b,sum,exact;

      a=0.0;
      b=1.0;
      n=20;
      exact=1.0/3.0;
      composite_midpoint(n,fx_sample,a,b,sum);
      std::cout << "approximation= " << sum  << '\n';
      std::cout << "exact= " << exact  << '\n';
      std::cout << "error= " << fabs(sum-exact)  << '\n';

      }

