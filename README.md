# ShearEmu

Source code for tangential shear emulator. The emulator has been trained
to predict the tangential shear from an N-body simulation based on on a 
WMAP7 cosmology. It uses Gaussian Process regression to calculate the fully non-linear 
P_gm(k) and xi_gm(r), with the galaxy clustering modelled by a 5 parameter Halo Occupation 
Distribution Model. The expressions for P_gm(k) and xi_gm(r) can then be used to 
derive the tangential shear using gamma_t.c. 

Please modify the supplied Makefile for your setup. You will need to 
link to GSL (any version with spline and linear algebra capabilities is fine).

To run the emulator, you will need to supply a parameter file containing the HOD
parameters and the desired output redshift. 
Please see the included example "params.ini" for the format. You will also
need to pass in the name of an output file on the command line.

Please reference this paper: https://adsabs.net/abs/2015ApJ...810...35K
if you use my code. 
