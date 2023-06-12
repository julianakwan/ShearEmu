#include<math.h>
#include <stdio.h>
#include <stdlib.h>

double func(double x, double a0, double b1, double b2, double b3); 
double func2(double x, double a0, double a1, double a2);




int pade(int nk, double *k, double *Pk, int npts, double *k_out, double *Pk_out)
{
// This code calculates the pade approximant which is used to extrapolate the power spectrum
// to higher k (beyond the emulator) for conversion to the two-point correlation function in configuration space
//
// Linear bias is assumed on large scales (k > 0.001), then the emulator is used and then the pade approximant
  
  int i; 
  double x0, y0, x1, y1, x2, y2, x3, y3; 


  x0 = k[nk-50]; y0 = Pk[nk-50];
  x1 = k[nk-25]; y1 = Pk[nk-25];
  x2 = k[nk-1];  y2 = Pk[nk-1];



  double a0 = 0;
  double b1 = 0;
  double b2 = 0;
  double b3 = 0; 

  double a1 = 0; 
  double a2 = 0; 

  /// Approximation of the form: y = a0*x*x+a1*x+a2
  a1  = y1 - y2 + (y1-y0)*(x1*x1-x2*x2)/(x0*x0-x1*x1);
  a1 /= ((x1-x0)*(x1*x1-x2*x2)/(x0*x0-x1*x1)+x1-x2);
  
/*   a0  = y0-y1+a1*(x1-x0); */
/*   a0 /= (x0*x0-x1*x1); */

  a2 = y0 - a0*x0*x0 - a1*x0;

  ///// Approximation of the form: y = a0/(1+b1*x+b2*x*x)//////

/*   b2 = (y1 - y0) + (y1 - y2)*(y0*x0-y1*x1)/(y1*x1-y2*x2); */
/*   b2 /= ((y0*x0*x0-y1*x1*x1)+(y2*x2*x2-y1*x1*x1)*(y0*x0-y1*x1)/(y1*x1-y2*x2)); */

/*   b1 = (y2-y1)+b2*(y2*x2*x2-y1*x1*x1); */
/*   b1 /= (y1*x1-y2*x2); */

/*   a0 = y0*(1.+b1*x0+b2*x0*x0); */

  ////// Approximation of the form:  y = a0/(1+b2*x*x)  ///////
/*   b2 = (y0 - y1)/(x1*x1*y1-x0*x0*y0); */
/*   a0 = y0*(1.+b2*x1*x1); */


  //// Approximation of the form: y = a0/(1+b1*x+b2*x*x+b3*x*x*x);
/*   double reslt = x1*x1*y1-x3*x3*y3+(x2*x2*y2-x3*x3*y3)*(x3*y3-x1*y1)/(x2*y2-x3*y3); */

/*   b3 = reslt*(y0-y1)+y3-y1+(x3*y3-x1*y1)*(y3-y2)/(x2*y2-x3*y3); */
/*   b3 /= (reslt*(x1*x1*x1*y1-x0*x0*x0*y0+(x1*y1-x0*y0)/(x2*y2-x3*y3)*(x3*x3*x3*y3-x2*x2*x2*y2))+(x3*y3-x1*y1)/(x2*y2-x3*y3)*(x2*x2*x2*y2-x3*x3*x3*y3)+x1*x1*x1*y1-x3*x3*x3*y3); */

/*   b2 = y3-y1+(x3*y3-x1*y1)/(x2*y2-x3-y3)*(y3-y2+b3*(x3*x3*x3*y3-x2*x2*x2*y2))+b3*(x3*x3*x3*y3-x1*x1*x1*y1); */
/*   b2 /= reslt; */

/*   b1 = y3-y2 + b2*(x3*x3*y3-x2*x2*y2)+b3*(x3*x3*x3*y3-x2*x2*x2*y2); */
/*   b1 /= x2*y2-x3*y3; */

/*   a0 = y0+b1*x0*y0+b2*x0*x0*y0+b3*x0*x0*x0*y0; */



  //Straight up linear bias on large scales
  for (i = 0; i < 379; i++)
    {
      Pk_out[i] = Pk[0]; 
      //      Pk_out[i] = 1.;
    }
  // Leave emulated regime untouched
  for(i=0; i< 330; i++)
    {
      Pk_out[i+379]=Pk[i]; 
      //      Pk_out[i] = 1.;
    }
  // Now do Pade approximation
  for(i = 709; i < npts; i++)
    {
      double k = pow(10.,k_out[i]);
      //      double k = k_out[i];
      //      double bias_ext = func(k, a0, b1, b2, b3); 
      double bias_ext= func2(k, a0, a1, a2); 
      Pk_out[i] = bias_ext*exp(-k/50); 
      //      Pk_out[i] = 1.;
    }


  return(0);
}

double func(double x, double a0, double b1, double b2, double b3)
{

  double f = a0/(1+b1*x+b2*x*x); 

  return(f); 
}

double func2(double x, double a0, double a1, double a2)
{

  double f = a0*x*x+a1*x+a2;

  return(f); 
}

