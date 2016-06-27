/*
  
 */
#define _USE_MATH_DEFINES //used for the value of pi in double precision(M_PI)
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//using namespace std;


void func_to_graph(double x,double& f) {
 
  f=sin(x);
  return;

}

int main() {
	
  double a,b,pi,h,x,func, func2;
  int N,iplot;
  pi=4.0*atan(1.0); //value of pi
  N=100; //number of points to plot
  a=0.0; //
  b=2.0*pi; //2pi
  ofstream graphfile,
           graphfile2;
  graphfile.open("file_to_graph");
  graphfile2.open("file_to_graph2");
  h=(b-a)/N;
  for (iplot=0;iplot<=N;iplot++) {
    x=a+h*iplot;
    func_to_graph(x,func);     
    if (x < pi) func2 = -1;
    else if (x >= pi) func2 = 1;
    graphfile << x << ' ' << func << '\n'; //func stores the value of y=f(x)
    graphfile2 << x << ' ' << func2 << '\n';
  }
  graphfile.close();     
  graphfile2.close();
  std::cout << "The file to graph is named: `file_to_graph' \n";
  return 0;
}

