#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h> //<cstdio> makes you  use std::
#include <stdlib.h>

using namespace std;

      void nfactorial(int n,double& nfact) {
   
      int i;

       nfact=1.0;
       for (i=1;i<=n;i++) {
        nfact=nfact*i;
       }

      }

      double power(double x,int N) {
 
      double pow=1.0;
      for (int i=1;i<=N;i++)
       pow*=x;

      return pow;

      }

      int main() {

      int N;
      double PN,exact,term,x,nfact;

      N=5;
      x=0.1;
      PN=1.0;
      for (int i=1;i<=N;i++) {
       nfactorial(i,nfact);
       term=power(x,i)/nfact;
       PN=PN+term;
      } 
      exact=exp(x);
      std::cout << "x= " << x << " and N= " << N << '\n';
      std::cout << "Taylor series approximation is: " << PN << '\n';
      std::cout << "exp(x) is: " << exact << '\n';
      std::cout << "error is: " << fabs(exact-PN) << '\n';

      }

