# ShearEmu

Source code for tangential shear emulator. Calculates P_gm(k) and
xi_gm(r) using a 5 parameter Halo Occupation Distribution Model to
model the galaxy catalog - these can then be used as inputs for
g-g lensing. 

Please use the supplied Makefile. You will need GSL (any version with
spline and linear algebra capabilities is fine).

To run the emulator, you will need to supply a parameter file. Please
see the included example "params.ini" for the format. You will also
need to pass in the name of an output file on the command line.


