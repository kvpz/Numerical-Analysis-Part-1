#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
      
using namespace std;

       // y=x^N  ln y=N ln x   y=e^{N ln x}
      double realpower(double x,double N) {
 
       return exp(N*log(x));

      }

      void fsample(double t,double& ff) {

       ff=realpower(t,5.0)-2.0;

      }

      void fsampleprime(double t,double& ffprime) {

       ffprime=5.0*realpower(t,4.0);

      }


      int main() { 

      double tolerance;
      int errormet,numiter;
      double t_old,f_old,f_oldprime,t_new,f_new;

      tolerance=1.0e-10;

      numiter=0;
      errormet=0;
      t_old=1.0;
      while (errormet==0) {

        fsample(t_old,f_old);
        fsampleprime(t_old,f_oldprime);
        t_new=t_old-f_old/f_oldprime;

        if (fabs(t_new-t_old)<tolerance) 
          errormet=1;
        t_old=t_new;
        numiter++;
        fsample(t_new,f_new); 
        std::cout << "numiter, t_new, f(t_new) " << numiter << ' ' << 
         t_new << ' ' << f_new << '\n'; 

      }
      fsample(t_new,f_new); 
      std::cout << "numiter, t_new, f(t_new) " << numiter << ' ' << 
         t_new << ' ' << f_new << '\n'; 

      }
