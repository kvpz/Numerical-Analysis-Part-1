#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// project 5
// implement piecewise cubic Hermite polynomial interpolation
// test on "four-leaf clover shape"
// show graphs, check error for increasing n.  Compare with spline interp.
//
// extra credit: implement a "periodic tridiagonal solve" instead of the
// inefficient partial pivoting solver used now.

   
        void matrix_solve(double** AA,double* xx,double* bb,int& status,
         int numelem) {
        
        double alpha,holdvalue;
        int i,j,k,holdj;

        status=1;
        for (i=1;i<=numelem-1;i++) {
         holdj=i;
         holdvalue=fabs(AA[i-1][i-1]);
         for (j=i+1;j<=numelem;j++) {
          if (fabs(AA[j-1][i-1])>holdvalue) {
           holdj=j;
           holdvalue=fabs(AA[j-1][i-1]);
          }
         }
         if (holdj!=i) {
          for (j=i;j<=numelem;j++) {
           holdvalue=AA[i-1][j-1];
           AA[i-1][j-1]=AA[holdj-1][j-1];
           AA[holdj-1][j-1]=holdvalue;
          }
         }
         holdvalue=bb[i-1];
         bb[i-1]=bb[holdj-1];
         bb[holdj-1]=holdvalue;
         if (fabs(AA[i-1][i-1])<1.0E-32) 
          status=0;
         else {
          for (j=i+1;j<=numelem;j++) {
           alpha=AA[j-1][i-1]/AA[i-1][i-1];
           for (k=i;k<=numelem;k++)
            AA[j-1][k-1]=AA[j-1][k-1]-alpha*AA[i-1][k-1];
           bb[j-1]=bb[j-1]-alpha*bb[i-1];
          }
         }
        }

        for (i=numelem;i>=1;i--) {
         if (status!=0) {
          holdvalue=bb[i-1];
          for (j=i+1;j<=numelem;j++)
           holdvalue=holdvalue-AA[i-1][j-1]*xx[j-1];
          if (fabs(AA[i-1][i-1])<1.0E-32) 
           status=0;
          else
           xx[i-1]=holdvalue/AA[i-1][i-1];
         }
        }

        }

      double xt_twoleaf(double t) {

      double tper;

      tper=t;
      while (tper<0.0) {
       tper=tper+4.0;
      }
      while (tper>4.0) {
       tper=tper-4.0;
      }
      double pi_local=4.0*atan(1.0);
      return (1.0+0.5*cos(pi_local*tper))*cos(pi_local*tper/2.0);

      }

      double yt_twoleaf(double t) {

      double tper;

      tper=t;
      while (tper<0.0) {
       tper=tper+4.0;
      }
      while (tper>4.0) {
       tper=tper-4.0;
      }
      double pi_local=4.0*atan(1.0);
      return (1.0+0.5*cos(pi_local*tper))*sin(pi_local*tper/2.0);

      }



      double xt_square(double t) {

      double tper;

      tper=t;
      while (tper<0.0) {
       tper=tper+4.0;
      }
      while (tper>4.0) {
       tper=tper-4.0;
      }

      if (tper<=1.0) 
       return tper;
      else if (tper<=2.0)
       return 1.0;
      else if (tper<=3.0) 
       return 1.0-(tper-2.0);
      else
       return 0.0;

      }

      double yt_square(double t) {

      double tper;

      tper=t;
      while (tper<0.0) {
       tper=tper+4.0;
      }
      while (tper>4.0) {
       tper=tper-4.0;
      }

      if (tper<=1.0) 
       return 0.0;
      else if (tper<=2.0)
       return tper-1.0;
      else if (tper<=3.0) 
       return 1.0;
      else
       return 1.0-(tper-3.0);

      }


      void get_pl_coefs(double* xx,double* ff,int n,
       double* acoef,double* bcoef) {

      int i;

      for (i=0;i<=n-1;i++) {
       acoef[i]=ff[i];
       bcoef[i]=(ff[i+1]-ff[i])/(xx[i+1]-xx[i]);
      }

      }


// this is very ineffecient - for many extra credit points
// write a much more efficient algorithm for finding the periodic
// spline coefficients. (a tridiagonal solve, but for the periodic
// case that arises here)
//
      void get_spline_coefs(double* xx,double* ff,int n, 
       double* acoef,double* bcoef,double* ccoef,double* dcoef) {

      double* harray=new double[n];

      int i,j,status;
      double** AA=new double*[n];
      for (i=0;i<n;i++)
       AA[i]=new double[n];
      double* bb=new double[n];
      double* cc=new double[n];

      for (i=0;i<=n;i++)
       acoef[i]=ff[i];

      for (i=0;i<n;i++)
       harray[i]=xx[i+1]-xx[i];
 
      for (i=0;i<n;i++)
      for (j=0;j<n;j++)
       AA[i][j]=0.0;
 
      for (i=1;i<=n;i++) {
       if (i==1) {
        AA[i-1][i-1]=2.0*(harray[i-1]+harray[n-1]);
        AA[i-1][i]=harray[i-1];
        AA[i-1][n-1]=harray[n-1];
        bb[i-1]=3.0*( (ff[1]-ff[0])/harray[0]- 
         (ff[n]-ff[n-1])/harray[n-1] );
       } else if (i==n) {
        AA[i-1][i-1]=2.0*(harray[n-1]+harray[n-2]);
        AA[i-1][i-2]=harray[n-2];
        AA[i-1][0]=harray[n-1];
        bb[i-1]=3.0*( (ff[i]-ff[i-1])/harray[i-1]- 
         (ff[i-1]-ff[i-2])/harray[i-2] );
       } else {
        AA[i-1][i-1]=2.0*(harray[i-1]+harray[i-2]);
        AA[i-1][i-2]=harray[i-2];
        AA[i-1][i]=harray[i-1];
        bb[i-1]=3.0*( (ff[i]-ff[i-1])/harray[i-1]- 
         (ff[i-1]-ff[i-2])/harray[i-2] );
       }
      }

      matrix_solve(AA,cc,bb,status,n);
      if (status!=1)
       std::cout << "matrix solve failed\n"; 
      for (i=0;i<n;i++)
       ccoef[i]=cc[i];
      ccoef[n]=ccoef[0];
 
      for (i=0;i<n;i++) {
       bcoef[i]=(acoef[i+1]-acoef[i])/harray[i]- 
          harray[i]*(ccoef[i+1]+2.0*ccoef[i])/3.0;
       dcoef[i]=(ccoef[i+1]-ccoef[i])/(3.0*harray[i]);
      }
    
      for (i=0;i<n;i++)
       delete AA[i];
      delete AA;
      delete bb;
      delete cc;

      } 


      void eval_pl_interp(double* xx,int n,double* acoef,
       double* bcoef,double x,double& plval) {

      int i;

      if (x<xx[0]) 
       plval=acoef[0];
      else if (x>xx[n]) 
       plval=acoef[n-1]+(xx[n]-xx[n-1])*bcoef[n-1];
      else {
       i=0;
       while ((x>xx[i+1])&&(i<n-1))
        i++;
       plval=acoef[i]+bcoef[i]*(x-xx[i]);
      }

      }
   

      void eval_spline_interp(double* xx,int n,double* acoef,
       double* bcoef,double* ccoef,double* dcoef,double x,
       double& splineval) {

      double dx;
      int i;

      if (x<xx[0]) 
       splineval=acoef[0];
      else if (x>xx[n]) {
       dx=xx[n]-xx[n-1];
       splineval=acoef[n-1]+bcoef[n-1]*dx+ 
        ccoef[n-1]*(dx*dx)+dcoef[n-1]*(dx*dx*dx);
      } else {
       i=0;
       while ((x>xx[i+1])&&(i<n-1)) 
        i++;
       dx=x-xx[i];
       splineval=acoef[i]+bcoef[i]*dx+ccoef[i]*(dx*dx)+ 
         dcoef[i]*(dx*dx*dx);
      }

      }

      int main() {

      int n=20;  // n is the number of intervals; n+1 is the number of nodes.
      int interp_type=2;  // 1=piecewise linear, 2=periodic cubic spline 

        // these are function variables (really neat feature of c++)
      double (*xt)(double);
      double (*yt)(double);

      double tlo=0.0;
      double thi=4.0;

    xt=&xt_square;
    yt=&yt_square;
//      xt=&xt_twoleaf;
//      yt=&yt_twoleaf;

      double* tsten=new double[n+1];
      double* xsten=new double[n+1];
      double* ysten=new double[n+1];

      double* acoefx=new double[n+1];
      double* bcoefx=new double[n+1];
      double* ccoefx=new double[n+1];
      double* dcoefx=new double[n+1];

      double* acoefy=new double[n+1];
      double* bcoefy=new double[n+1];
      double* ccoefy=new double[n+1];
      double* dcoefy=new double[n+1];

      double h,hplot;
      double tplot,xplot,yplot;
      double xexact,yexact,L2_error;
      int i,nplot;

      h=(thi-tlo)/n;
      for (i=0;i<=n;i++) {
       tsten[i]=tlo+i*h;
       xsten[i]=xt(tsten[i]);
       ysten[i]=yt(tsten[i]);
      }

      if (interp_type==1) {
       get_pl_coefs(tsten,xsten,n,acoefx,bcoefx);
       get_pl_coefs(tsten,ysten,n,acoefy,bcoefy);
      } else if (interp_type==2) {
       get_spline_coefs(tsten,xsten,n, 
          acoefx,bcoefx,ccoefx,dcoefx);
       get_spline_coefs(tsten,ysten,n, 
          acoefy,bcoefy,ccoefy,dcoefy);
      } else
       std::cout <<"interp_type invalid\n";

      ofstream interpfile;
      interpfile.open("curveinterp");
      std::cout << "output filename is curveinterp \n";
      std::cout << "plot in gnuplot with option: using 2:3\n";

      nplot=1000;
      hplot=(thi-tlo)/nplot;

      std::cout << "tlo,thi,n,h " << tlo << ' ' << thi << ' ' <<
         n << ' ' << h << '\n';
      std::cout << "tlo,thi,nplot,hplot " << tlo << ' ' << thi << ' ' <<
       nplot << ' ' << hplot << '\n';

      L2_error=0.0;

      for (i=0;i<=nplot;i++) {

       tplot=tlo+i*hplot;
  
       if (interp_type==1) {
        eval_pl_interp(tsten,n,acoefx,bcoefx,tplot,xplot);
        eval_pl_interp(tsten,n,acoefy,bcoefy,tplot,yplot);
       } else if (interp_type==2) {
        eval_spline_interp(tsten,n, 
          acoefx,bcoefx,ccoefx,dcoefx,tplot,xplot);
        eval_spline_interp(tsten,n, 
          acoefy,bcoefy,ccoefy,dcoefy,tplot,yplot);
       } else
        std::cout << "interp_type invalid\n";

       xexact=xt(tplot);
       yexact=yt(tplot);
       L2_error=L2_error+((xexact-xplot)*(xexact-xplot)+
                          (yexact-yplot)*(yexact-yplot))*hplot;
       interpfile << tplot << ' ' << xplot << ' ' << yplot << '\n';        

      }
      interpfile.close();

      L2_error=sqrt(L2_error);

      if (interp_type==1) 
       std::cout << "interp_type=1 (piecewise linear interpolation)\n";
      else if (interp_type==2) 
       std::cout << "interp_type=2 (periodic cubic spline interpolation)\n";
      else
       std::cout << "interp_type out of range\n";

      std::cout << "n and h are " << n << ' ' << h << '\n';
      std::cout << "L2 error is " << L2_error << '\n';

      delete tsten;
      delete xsten;
      delete ysten;
      delete acoefx;
      delete bcoefx;
      delete ccoefx;
      delete dcoefx;
      delete acoefy;
      delete bcoefy;
      delete ccoefy;
      delete dcoefy;

      }
  
