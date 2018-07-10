#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"gsl/gsl_spline.h"
#include"gsl/gsl_integration.h"
#include"gsl/gsl_sf_bessel.h"
#include"gamma_t.h"
#include "k_m_ext.h"
#include "params.h"


int read_nz_dist(char filename[], int ncol, int *n, double *z, double *chi, double *n_z); 
int read_nz_dist_2(char filename[], int ncol, int *n, double *z, double *chi, double *n_z); 
double integrate_source_dist(double chi, int nbins, double *z, double *chi_dash, double *n_z); 
double comoving_int(double z, void *params); 
double angular_diameter_distance(struct cosmo_params cosmo, double redshift); 
int z_to_comoving(int nbins, double *z, double *chi, struct cosmo_params cosmo); 
double source_dist_integrand (double chi_dash, void * params); 
int setup_hankel_transform(struct info_for_hankel *hankel); 
int setup_gamma_t(int nchi, int ntheta, double *theta, struct survey_info survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info, struct cosmo_params cosmo); 
int setup_DeltaSigma(int nchi, int ntheta, double *theta, struct survey_info survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info,  struct cosmo_params cosmo); 
int clean_up_gamma_t(struct survey_info *survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info); 
int clean_up_DeltaSigma(struct survey_info *survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info);
int hankel_transform_gamma_t(struct info_for_hankel hankel, struct cosmo_params cosmo, int nk, double *k, double *Pk, int nr, double *theta, double chi, double *out); 
int hankel_transform_DeltaSigma(struct info_for_hankel hankel, struct cosmo_params cosmo, int nk, double *k, double *Pk, int nr, double *theta, double chi, double *out); 
int integrate_lens_dist_gamma_t(struct survey_info survey, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel, struct info_for_lens_dist_integrand info, 
			int ntheta, double *gamma_t); 
int integrate_lens_dist_DeltaSigma(struct survey_info survey, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel, struct info_for_lens_dist_integrand info, 
			int ntheta, double *gamma_t); 
double lens_dist_integrand(double chi, void *params); 
int tabulate_fk_integral(int nbins_l, double *z_l, double *chi_l, double *n_l, int nbins_s, double *z_s, double *chi_s, double *n_s, int nchi, double *table); 
int tabulate_pk_integral_gamma_t(struct survey_info survey, struct info_for_lens_dist_integrand info, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel); 
int tabulate_pk_integral_DeltaSigma(struct survey_info survey, struct info_for_lens_dist_integrand info, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel); 
int normalise_distributions(int n, double *z, double *chi, double *n_z, struct cosmo_params cosmo); 
int convert_params(struct HOD_params HOD, double params[]); 
int calc_gamma_t(char source_filename[], char lens_filename[], struct HOD_params HOD, struct cosmo_params cosmo); 
int calc_DeltaSigma(char source_filename[], char lens_filename[], struct HOD_params HOD, struct cosmo_params cosmo); 

// Constants
double c = 2997.92458;///0.71; // Mpc/h
double c2_over_G = 2.092e19; // M_sun/Mpc
double arcminutes_to_radians = M_PI/180./60.;

int main(int argc, char**argv)
{
  struct cosmo_params cosmo; 

  cosmo.Om = 0.2648;
  cosmo.Omb = 0.04479;
  cosmo.w = -1;
  cosmo.ns = 0.963;
  cosmo.h = 0.71;
  cosmo.sigma8 = 0.8;
  cosmo.redshift = 0.;
  cosmo.growth_factor = 1.;

//  HOD parameters
  struct HOD_params HOD; 


  HOD.Mcut = 13.; 
  HOD.M1 = 14.;
  HOD.sigma = 1.0; 
  HOD.kappa = 1.0; 
  HOD.alpha = 1.0;




  char source_filename[256], lens_filename[256];

  // Read in source distribution
  sprintf(source_filename, "/users/astro/jkwan/ggl/data/id_im3shape_z0.50-1.30.txt_nofz_tpz");
  //  sprintf(source_filename, "/home/kjuliana/ggl/ggl_source_redshift_distribution.dat");
  /* sprintf(source_filename, "/Users/julianakwan/Documents/DES/ggl/id_im3shape_z0.50-1.30.txt_nofz_tpz"); */


  // Read in lens distribution
   sprintf(lens_filename, "/users/astro/jkwan/ggl/data/id_redmagic_z0.20-0.40_L1.00.txt_nofz_tpz");
  /* sprintf(lens_filename, "/home/kjuliana/ggl/ggl_lens_redshift_distribution.dat"); */
  /* sprintf(lens_filename, "/Users/julianakwan/Documents/DES/ggl/id_redmagic_z0.20-0.40_L1.00.txt_nofz_tpz"); */


  calc_gamma_t(source_filename, lens_filename, HOD, cosmo);  

  calc_DeltaSigma(source_filename, lens_filename, HOD, cosmo);  

  return(0);
}


int read_nz_dist(char filename[], int ncol, int *n, double *z, double *chi, double *n_z)
{
  FILE *fp = fopen(filename, "r"); 
  int nread = ncol; 

  float any[2]; 

  (*n) = 0; 
  while (nread == ncol)
    {
      //      nread = fscanf(fp, "%lf %lf %lf %lf", &(z[*n]), &(chi[*n]), &(any[1]), &(n_z[*n])); 
      nread = fscanf(fp, "%lf %lf", &(z[*n]), &(n_z[*n])); 
      (*n)++; 
    }
  (*n)--; 

  
  fclose(fp); 
  return(0); 
}

int read_nz_dist_2(char filename[], int ncol, int *n, double *z, double *chi, double *n_z)
{
  FILE *fp = fopen(filename, "r"); 
  int nread = ncol; 

  float any[2]; 

  (*n) = 0; 
  while (nread == ncol)
    {
      nread = fscanf(fp, "%lf %lf", &(z[*n]), &(n_z[*n])); 
      (*n)++; 
    }
  (*n)--; 

  //normalise distributions
  /* gsl_interp_accel *acc  = gsl_interp_accel_alloc (); */
  /* gsl_spline *spline     = gsl_spline_alloc (gsl_interp_cspline, *n); */
  /* gsl_spline_init (spline, z, n_z, *n); */

  /* double norm = gsl_spline_eval_integ(spline, z[0], z[*n-1], acc);  */

  /* gsl_spline_free (spline); */
  /* gsl_interp_accel_free (acc); */

  /* //  fprintf(stderr, "%lf\n", norm);  */

  /* int i;  */

  /* for (i=0; i < (*n); i++) */
  /*   { */
  /*     n_z[i]/=norm; */
  /*     //      fprintf(stderr, "%lf %7.5e\n", z[i], n_z[i]); */
  /*   } */
  
  fclose(fp); 
  return(0); 
}

int normalise_distributions(int n, double *z, double *chi, double *n_z, struct cosmo_params cosmo)
{
  int i;
  
  for (i = 0; i < n; i++)
    {
      double Hz = (cosmo.Om*(1+z[i])*(1+z[i])*(1+z[i])+(1-cosmo.Om)); 
      n_z[i] /= sqrt(Hz); 
      /* n_z[i] *= c; */
    }

  gsl_interp_accel *acc  = gsl_interp_accel_alloc ();
  gsl_spline *spline     = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline_init (spline, chi, n_z, n);

  double norm = gsl_spline_eval_integ(spline, chi[0], chi[n-1], acc);

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  fprintf(stderr, "%lf\n", norm);



  for (i=0; i < n; i++)
    {
      n_z[i]/=norm;
      //   fprintf(stderr, "%lf %7.5e\n", chi[i], n_z[i]);
    }

  return(0); 
}

double integrate_source_dist(double chi, int nbins, double *z, double *chi_dash, double *n_z)
{
  // this part does this integral f(\chi') = \int d\chi' D_l*D_ls/Ds/a*n_s(\chi')

  double f_chi; 

  // spline chi(z) to get redshifts at chi, chi'
  gsl_interp_accel *acc_chi = gsl_interp_accel_alloc(); 
  gsl_spline *spline_chi = gsl_spline_alloc(gsl_interp_cspline, nbins); 

  gsl_spline_init(spline_chi, chi_dash, z, nbins); 

  // struct of parameters that need to be passed to the 
  // function to evaluate the integrand
  struct info_for_source_dist_integrand info; 

  info.integrand = malloc(3*nbins*sizeof(double)); 

  int i; 

  for (i=0; i<nbins; i++)
    {
      info.integrand[3*i+0] = chi_dash[i]; 
      info.integrand[3*i+1] = z[i]; 
      info.integrand[3*i+2] = n_z[i]; 
    }
  info.npts = nbins; 
  info.chi = chi; // comoving coordinate
  //  info.z_l = gsl_spline_eval(spline_chi, chi, acc_chi);  

  double error;

  gsl_integration_cquad_workspace *w 
    = gsl_integration_cquad_workspace_alloc(1000);

  
  gsl_function F;
  F.function = &source_dist_integrand;
  F.params = &info;

  // check that integrand is non-zero
  // find last bin in n_z_source that contains non-zero entry



  if (chi > chi_dash[nbins-1])
    {
      f_chi = 0;
    }
  else
    {
      gsl_integration_cquad(&F, chi, chi_dash[nbins-1], 0, 1e-5, w, &f_chi, &error, 0);
      //      gsl_integration_qags (&F, chi, chi_dash[nbins-1], 0, 1e-5, 10000,  w, &(f_chi), &error); 
      //      fprintf(stderr, "HERE: %lf %lf\n", chi, f_chi); 
    }
      //      fprintf(stderr, "Here\n"); 
  //  f_chi *= chi; 

  gsl_integration_cquad_workspace_free(w);

  gsl_spline_free (spline_chi);
  gsl_interp_accel_free (acc_chi);
  free(info.integrand); 

  return(f_chi); 
}


double source_dist_integrand (double chi_dash, void * params) 
{
  // This is D_l * D_ls / D_s / a(z_l) * W_s(chi)
  struct info_for_source_dist_integrand info = *(struct info_for_source_dist_integrand*) (params); 

  int i = 0, npts=info.npts; 
  double *x = malloc(npts*sizeof(double));  // chi
  double *y = malloc(npts*sizeof(double));  // n_z (sources)
  /* double *z = malloc(npts*sizeof(double));  //redshifts */

  for (i = 0; i < npts; i++)
    {
      x[i] = info.integrand[3*i+0]; 
      y[i] = info.integrand[3*i+2]; 
      /* z[i] = info.integrand[3*i+1];  */
    }
  // Need to pass in chi, z_l, npts  


  gsl_interp_accel *acc = gsl_interp_accel_alloc(); 
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, npts); 

  gsl_spline_init(spline, x, y, npts); 

  /* gsl_interp_accel *acc_z = gsl_interp_accel_alloc();  */
  /* gsl_spline *spline_z = gsl_spline_alloc(gsl_interp_cspline, npts);  */

  /* gsl_spline_init(spline_z, x, z, npts);  */

  double n_s; 

   if (chi_dash >= x[0] && chi_dash <= x[npts-1])
       n_s = gsl_spline_eval(spline, chi_dash, acc);
   else
     n_s = 0;

  double chi = info.chi;
  /* double z_l = info.z_l;  */
  /* double z_s = gsl_spline_eval(spline_z, chi_dash, acc_z);  */
  /* double D_l = chi/(1.+z_l);  */
  /* double D_s = chi_dash/(1.+z_s);  */
  /* double D_ls = (chi_dash-chi)/(1.+(z_s-z_l));  */
  /* double a_z_l = 1./(1.+z_l);  */
  double f;

  if (chi_dash > chi)
    //    f = D_l*D_ls/D_s/a_z_l*n_s; 
    //    f = chi*(chi_dash-chi)/chi_dash*n_s;
    f = (chi_dash-chi)/chi_dash*n_s;
  else
    f = 0; 

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  free(x); free(y); 

  return f;
}

int calc_gamma_t(char source_filename[], char lens_filename[], struct HOD_params HOD, struct cosmo_params cosmo) 
{
  int i; 

  // Read in source and lens distributions

  struct survey_info survey; 
  survey.nbins_s = 0;
  survey.z_s = malloc(1000*sizeof(double));
  survey.n_z_s = malloc(1000*sizeof(double));
  survey.chi_s = malloc(1000*sizeof(double));

  read_nz_dist(source_filename, 2, &survey.nbins_s, survey.z_s, survey.chi_s, survey.n_z_s);
  z_to_comoving(survey.nbins_s, survey.z_s, survey.chi_s, cosmo);
  normalise_distributions(survey.nbins_s, survey.z_s, survey.chi_s, survey.n_z_s, cosmo); 

  survey.nbins_l = 0;
  survey.z_l = malloc(1000*sizeof(double));
  survey.n_z_l = malloc(1000*sizeof(double));
  survey.chi_l = malloc(1000*sizeof(double));


  read_nz_dist(lens_filename, 2, &survey.nbins_l, survey.z_l, survey.chi_l, survey.n_z_l);
  z_to_comoving(survey.nbins_l, survey.z_l, survey.chi_l, cosmo);
  normalise_distributions(survey.nbins_l,  survey.z_l, survey.chi_l, survey.n_z_l, cosmo); 


  // Set up theta values
  int ntheta = 50; 
  double theta_min = 0.1;
  double theta_max = 100.;

  double *theta = malloc(ntheta*sizeof(double));
  double *model = malloc(ntheta*sizeof(double));

  for (i = 0; i < ntheta; i++)
    {
      theta[i] = log10(theta_min) + (double)i*(log10(theta_max)-log10(theta_min))/((double)ntheta-1); 
      theta[i] = pow(10.,theta[i]); 
    }


  static struct info_for_hankel hankel;
  static struct info_for_lens_dist_integrand info;
  int n_emu_redshifts = 6;


  setup_gamma_t(n_emu_redshifts, ntheta, theta, survey, &hankel, &info, cosmo);
  integrate_lens_dist_gamma_t(survey, HOD, cosmo, hankel, info, ntheta, model);

/*   setup_DeltaSigma(n_emu_redshifts, ntheta, theta, survey, &hankel, &info,cosmo); */
/*   integrate_lens_dist_DeltaSigma(survey, HOD, cosmo, hankel, info, ntheta, model); */
  
  // output values
  FILE *fp_check = fopen("gamma_t_emu_test.dat", "w");
  
  for (i = 0; i < ntheta; i++)
      fprintf(fp_check, "%e %e\n", theta[i], model[i]);

  fclose(fp_check);

  free(theta); free(model);

  // clean up 
  clean_up_gamma_t(&survey, &hankel, &info); 

  return(0); 
}

int calc_DeltaSigma(char source_filename[], char lens_filename[], struct HOD_params HOD, struct cosmo_params cosmo) 
{
  int i; 

  // Read in source and lens distributions

  struct survey_info survey; 
  survey.nbins_s = 0;
  survey.z_s = malloc(1000*sizeof(double));
  survey.n_z_s = malloc(1000*sizeof(double));
  survey.chi_s = malloc(1000*sizeof(double));

  read_nz_dist(source_filename, 2, &survey.nbins_s, survey.z_s, survey.chi_s, survey.n_z_s);
  z_to_comoving(survey.nbins_s, survey.z_s, survey.chi_s, cosmo);
  normalise_distributions(survey.nbins_s, survey.z_s, survey.chi_s, survey.n_z_s, cosmo); 

  survey.nbins_l = 0;
  survey.z_l = malloc(1000*sizeof(double));
  survey.n_z_l = malloc(1000*sizeof(double));
  survey.chi_l = malloc(1000*sizeof(double));


  read_nz_dist(lens_filename, 2, &survey.nbins_l, survey.z_l, survey.chi_l, survey.n_z_l);
  z_to_comoving(survey.nbins_l, survey.z_l, survey.chi_l, cosmo);
  normalise_distributions(survey.nbins_l,  survey.z_l, survey.chi_l, survey.n_z_l, cosmo); 


  // Set up theta values
  int nR = 50; 
  double R_min = 0.01;
  double R_max = 10.;

  double *R = malloc(nR*sizeof(double));
  double *DeltaSigma = malloc(nR*sizeof(double));

  for (i = 0; i < nR; i++)
    {
      R[i] = log10(R_min) + (double)i*(log10(R_max)-log10(R_min))/((double)nR-1); 
      R[i] = pow(10.,R[i]); 
    }


  static struct info_for_hankel hankel;
  static struct info_for_lens_dist_integrand info;
  int n_emu_redshifts = 6;


  setup_DeltaSigma(n_emu_redshifts, nR, R, survey, &hankel, &info, cosmo);
  integrate_lens_dist_DeltaSigma(survey, HOD, cosmo, hankel, info, nR, DeltaSigma);
  
  // output values
  FILE *fp_check = fopen("DeltaSigma_emu_test.dat", "w");
  
  for (i = 0; i < nR; i++)
      fprintf(fp_check, "%e %e\n", R[i], DeltaSigma[i]);

  fclose(fp_check);

  free(R); free(DeltaSigma);

  // clean up 
  clean_up_gamma_t(&survey, &hankel, &info); 

  return(0); 
}

int setup_gamma_t(int nchi, int ntheta, double *theta, struct survey_info survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info,  struct cosmo_params cosmo)
{
  int i, j; 

  (*info).nchi = nchi; // how many points in radial distance we want to precalculate.
  (*info).ntheta = ntheta; // how many points in theta we want

  (*info).theta = malloc(ntheta*sizeof(double)); 

  for (i = 0; i < ntheta; i++)
    (*info).theta[i] = theta[i]; 

  //////////// SET UP integral over source distribution (lens efficiency) /////////////////////////
  // set up table for f_K integral
  (*info).npts = survey.nbins_l; 
  (*info).chi = calloc((*info).npts,sizeof(double)); 


  for (i = 0; i < (*info).npts; i++)
    (*info).chi[i] = survey.chi_l[i]; 


  (*info).f_chi_table=malloc((*info).npts*sizeof(double)); 

   tabulate_fk_integral(survey.nbins_l, survey.z_l, survey.chi_l, survey.n_z_l, survey.nbins_s, survey.z_s, survey.chi_s,  survey.n_z_s, (*info).nchi, (*info).f_chi_table); 

  (*info).acc_f_K = gsl_interp_accel_alloc();
  (*info).spline_f_K = gsl_spline_alloc(gsl_interp_cspline, (*info).npts);

   gsl_spline_init((*info).spline_f_K, (*info).chi, (*info).f_chi_table, (*info).npts);

  /* (*info).n_l = malloc((*info).npts*sizeof(double));  */

  /* for (i = 0; i < survey.nbins_l; i++) */
  /*   (*info).n_l[i] = survey.n_z_l[i]/(1+survey.z_l[i])/(1+survey.z_l[i]);  */

  /* // set up spline for f_K integral */
  /* (*info).acc_f_K = gsl_interp_accel_alloc(); */
  /* (*info).spline_f_K = gsl_spline_alloc(gsl_interp_cspline, survey.nbins_l); */

  /* gsl_spline_init((*info).spline_f_K, (*info).chi, (*info).n_l, survey.nbins_l); // factor of (1+z) is there to change to proper coordinates */


  //////////// SET UP P(k) part of integral /////////////////////////
  // initialize values, J2, zeros etc. for Hankel transform 
  setup_hankel_transform(hankel); 

  // this is the chi that goes with the table of Pk values
  // they need to correspond to the redshifts of the emulator

  (*info).chi_table = calloc((*info).nchi, sizeof(double)); 


  double emu_scale_factors[6] = {0.9999,  0.9083, 0.8051, 0.6974,  0.6086, 0.4985}; 
  for (i = 0; i < (*info).nchi; i++)
    {
      /* info.chi_table[i] = log10(chi_min) + i*((log10(chi_max) - log10(chi_min))/((double)info.nchi-1)); */
      /* info.chi_table[i] = pow(10.,info.chi_table[i]); */

      /* info.chi_table[i] = chi_min + i*(chi_max-chi_min)/((double)info.nchi-1); */
      //      info.chi_table[i] /= c;

      //      (*info).chi_table[i] = survey.chi_l[0] + i*(survey.chi_l[survey.nbins_l-1]-survey.chi_l[0])/((double)(*info).nchi-1); 
      double emu_redshift = 1./emu_scale_factors[i] - 1; 
      (*info).chi_table[i] = angular_diameter_distance(cosmo, emu_redshift); 
    }

  (*info).pk_table = malloc((*info).nchi*sizeof(double*)); 
  for (i = 0; i < (*info).nchi; i++)
    (*info).pk_table[i] = calloc((*info).ntheta,sizeof(double)); 


  return(0); 
}

int setup_DeltaSigma(int nchi, int ntheta, double *theta, struct survey_info survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info,  struct cosmo_params cosmo)
{
  int i, j; 

  (*info).nchi = nchi; // how many points in radial distance we want to precalculate.
  (*info).ntheta = ntheta; // how many points in theta we want

  (*info).theta = malloc(ntheta*sizeof(double)); 

  for (i = 0; i < ntheta; i++)
    (*info).theta[i] = theta[i]; 

  //////////// SET UP integral over source distribution (lens efficiency) /////////////////////////
  // set up table for f_K integral
  (*info).npts = survey.nbins_l; 
  (*info).chi = calloc((*info).npts,sizeof(double)); 


  for (i = 0; i < (*info).npts; i++)
    (*info).chi[i] = survey.chi_l[i]; 


  (*info).f_chi_table=malloc((*info).npts*sizeof(double)); 

  (*info).n_l = malloc((*info).npts*sizeof(double)); 

  for (i = 0; i < survey.nbins_l; i++)
    (*info).n_l[i] = survey.n_z_l[i]/(1+survey.z_l[i])/(1+survey.z_l[i]); 

  // set up spline for f_K integral
  (*info).acc_f_K = gsl_interp_accel_alloc();
  (*info).spline_f_K = gsl_spline_alloc(gsl_interp_cspline, survey.nbins_l);

  gsl_spline_init((*info).spline_f_K, (*info).chi, (*info).n_l, survey.nbins_l); // factor of (1+z) is there to change to proper coordinates


  //////////// SET UP P(k) part of integral /////////////////////////
  // initialize values, J2, zeros etc. for Hankel transform 
  setup_hankel_transform(hankel); 

  // this is the chi that goes with the table of Pk values
  // they need to correspond to the redshifts of the emulator

  (*info).chi_table = calloc((*info).nchi, sizeof(double)); 


  double emu_scale_factors[6] = {0.9999,  0.9083, 0.8051, 0.6974,  0.6086, 0.4985}; 
  for (i = 0; i < (*info).nchi; i++)
    {
      /* info.chi_table[i] = log10(chi_min) + i*((log10(chi_max) - log10(chi_min))/((double)info.nchi-1)); */
      /* info.chi_table[i] = pow(10.,info.chi_table[i]); */

      /* info.chi_table[i] = chi_min + i*(chi_max-chi_min)/((double)info.nchi-1); */
      //      info.chi_table[i] /= c;

      //      (*info).chi_table[i] = survey.chi_l[0] + i*(survey.chi_l[survey.nbins_l-1]-survey.chi_l[0])/((double)(*info).nchi-1); 
      double emu_redshift = 1./emu_scale_factors[i] - 1; 
      (*info).chi_table[i] = angular_diameter_distance(cosmo, emu_redshift); 
    }

  (*info).pk_table = malloc((*info).nchi*sizeof(double*)); 
  for (i = 0; i < (*info).nchi; i++)
    (*info).pk_table[i] = calloc((*info).ntheta,sizeof(double)); 



  return(0); 
}


int clean_up_gamma_t(struct survey_info *survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info)
{
  int i; 

  free((*hankel).zero); free((*hankel).phi); free((*hankel).phidash); free((*hankel).weights); free((*hankel).besselterm); 

  for (i = 0; i < (*info).nchi; i++)
    free((*info).pk_table[i]);

  free((*info).pk_table);

  free((*info).f_chi_table); free((*info).chi); free((*info).chi_table); 

  gsl_spline_free((*info).spline_f_K); gsl_interp_accel_free((*info).acc_f_K); 

  free((*survey).z_l); free((*survey).n_z_l); free((*survey).chi_l); 
  free((*survey).z_s); free((*survey).n_z_s); free((*survey).chi_s); 


  return(0); 
}

int clean_up_DeltaSigma(struct survey_info *survey, struct info_for_hankel *hankel, struct info_for_lens_dist_integrand *info)
{
  int i; 

  free((*hankel).zero); free((*hankel).phi); free((*hankel).phidash); free((*hankel).weights); free((*hankel).besselterm); 

  for (i = 0; i < (*info).nchi; i++)
    free((*info).pk_table[i]);

  free((*info).pk_table); free((*info).chi_table); free((*info).chi); free((*info).n_l); 

  //  free((*info).f_chi_table); 

  gsl_spline_free((*info).spline_f_K); gsl_interp_accel_free((*info).acc_f_K); 

  free((*survey).z_l); free((*survey).n_z_l); free((*survey).chi_l); 
  free((*survey).z_s); free((*survey).n_z_s); free((*survey).chi_s); 

  return(0); 
}; 



int integrate_lens_dist_gamma_t(struct survey_info survey, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel, struct info_for_lens_dist_integrand info, int ntheta, double *gamma_t)
{
  // set up tables for P(k) integral
  int i,j; 

  tabulate_pk_integral_gamma_t(survey, info, HOD, cosmo, hankel);

  double error, ans;
  //  double *gamma_t = malloc(2*info.ntheta*sizeof(double)); 

  gsl_integration_cquad_workspace * w 
    = gsl_integration_cquad_workspace_alloc (100000);
  
  gsl_function F;

  /* double SigmaCrit = CalcSigmaCritAvg(survey, cosmo);  */
  /* fprintf(stderr, "%lf \n", SigmaCrit);  */
  for (j = 0; j < info.ntheta; j++)
    {

      // set up spline for pkJ integral
      info.acc_pkJ = gsl_interp_accel_alloc();
      info.spline_pkJ = gsl_spline_alloc(gsl_interp_cspline, info.nchi);


      double *pkJ = malloc(info.nchi*sizeof(double));
      for (i = 0; i < info.nchi; i++)
      	{
      	  pkJ[i] = info.pk_table[i][j];   // i = chi, j = theta
      	  //  fprintf(stderr, "%lf %e\n",info.chi_table[i], pkJ[i]);
      	}
      gsl_spline_init(info.spline_pkJ, info.chi_table, pkJ, info.nchi);


      F.function = &lens_dist_integrand;
      F.params = &info;

      //      gsl_integration_qag (&F, chi_l[0], chi_l[nbins_l-1], 0, 5e-4, 100000, GSL_INTEG_GAUSS61, w, &ans, &error);
      //      gsl_integration_qags (&F, chi_l[8], chi_l[20], 0, 1e-4, 100000,  w, &(ans), &error); // CHECK THIS chi_l[8] is the first non-zero entry

      //      gsl_integration_cquad(&F, survey.chi_l[0], survey.chi_l[survey.nbins_l-1], 0, 5e-5, w, &ans, &error, 0);
      gsl_integration_cquad(&F, survey.chi_l[0], survey.chi_l[91], 0, 5e-5, w, &ans, &error, 0); // survey.chi_l is zero after entry 91.

      /* gamma_t[2*j+0] = theta[j]; */
      /* gamma_t[2*j+1] = 3./2.*M_PI/c/c*cosmo.Om*ans; */

      gamma_t[j] = 3./2.*M_PI/c/c*cosmo.Om*ans;
      //      gamma_t[j] *= c2_over_G/1e6/1e6/cosmo.h/4./M_PI;
      //gsl_interp_accel_reset (info.acc_pkJ);
      free(pkJ);   
      gsl_spline_free(info.spline_pkJ); gsl_interp_accel_free(info.acc_pkJ);

    }


  gsl_integration_cquad_workspace_free (w);


  return(0); 
}

int integrate_lens_dist_DeltaSigma(struct survey_info survey, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel, struct info_for_lens_dist_integrand info, int ntheta, double *gamma_t)
{
  // set up tables for P(k) integral
  int i,j; 

  tabulate_pk_integral_DeltaSigma(survey, info, HOD, cosmo, hankel);

  double error, ans;
  //  double *gamma_t = malloc(2*info.ntheta*sizeof(double)); 

  gsl_integration_cquad_workspace * w 
    = gsl_integration_cquad_workspace_alloc (100000);
  
  gsl_function F;

  /* double SigmaCrit = CalcSigmaCritAvg(survey, cosmo);  */
  /* fprintf(stderr, "%lf \n", SigmaCrit);  */
  for (j = 0; j < info.ntheta; j++)
    {

      // set up spline for pkJ integral
      info.acc_pkJ = gsl_interp_accel_alloc();
      info.spline_pkJ = gsl_spline_alloc(gsl_interp_cspline, info.nchi);


      double *pkJ = malloc(info.nchi*sizeof(double));
      for (i = 0; i < info.nchi; i++)
      	{
      	  pkJ[i] = info.pk_table[i][j];   // i = chi, j = theta
      	  //  fprintf(stderr, "%lf %e\n",info.chi_table[i], pkJ[i]);
      	}
      gsl_spline_init(info.spline_pkJ, info.chi_table, pkJ, info.nchi);


      F.function = &lens_dist_integrand;
      F.params = &info;

      //      gsl_integration_qag (&F, chi_l[0], chi_l[nbins_l-1], 0, 5e-4, 100000, GSL_INTEG_GAUSS61, w, &ans, &error);
      //      gsl_integration_qags (&F, chi_l[8], chi_l[20], 0, 1e-4, 100000,  w, &(ans), &error); // CHECK THIS chi_l[8] is the first non-zero entry

      //      gsl_integration_cquad(&F, survey.chi_l[0], survey.chi_l[survey.nbins_l-1], 0, 5e-5, w, &ans, &error, 0);
      gsl_integration_cquad(&F, survey.chi_l[0], survey.chi_l[91], 0, 5e-5, w, &ans, &error, 0); // survey.chi_l is zero after entry 91.

      /* gamma_t[2*j+0] = theta[j]; */
      /* gamma_t[2*j+1] = 3./2.*M_PI/c/c*cosmo.Om*ans; */

      gamma_t[j] = 3./2.*M_PI/c/c*cosmo.Om*ans;
      gamma_t[j] *= c2_over_G/1e6/1e6/cosmo.h/4./M_PI;
      //gsl_interp_accel_reset (info.acc_pkJ);
      free(pkJ);   
      gsl_spline_free(info.spline_pkJ); gsl_interp_accel_free(info.acc_pkJ);

    }


  gsl_integration_cquad_workspace_free (w);


  return(0); 
}

int tabulate_pk_integral_gamma_t(struct survey_info survey, struct info_for_lens_dist_integrand info, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel)
{
  int i,j; 


  // set up arrays for emulator 

  int nk_emu = 330; 
  double *k_emu = malloc(nk_emu*sizeof(double)); 
  double *Pk_emu = malloc(nk_emu*sizeof(double)); 

  int extra_pts = 51; // no of points to extend in either direction 
  //  int nk_ext = nk_emu+2*extra_pts-2;
  int nk_ext = 976;
  double *k_ext  = malloc(nk_ext*sizeof(double));
  double *Pk_ext = malloc(nk_ext*sizeof(double));

  double params[5]; 
  convert_params(HOD, params); 

  for (i= 0; i < nk_emu; i++)
    {
      k_emu[i] = pow(10.,logk[i]); 
      k_emu[i]  /= cosmo.h;
      //      fprintf(stderr, "%lf\n", k_emu[i]); 
    }

  for (i = 0; i < nk_ext; i++)
    {
      k_ext[i] = pow(10.,k_m_ext[i]); 
      k_ext[i] /= cosmo.h;
    }

  double k_min = 1e-3;
  double k_max = 1e8;

  // These are determined by the snapshots used by the emulator: 
  double emu_scale_factors[6] = {0.9999,  0.9083, 0.8051, 0.6974,  0.6086, 0.4985};// last scalefactor is actually a = 1, but can't evaluate at chi = 0;  

  for (i = 0; i < info.nchi; i++)
    {
      // calculate the power spectrum from emulator

      cosmo.redshift = 1/emu_scale_factors[i] - 1; 

      emu(params, cosmo.redshift, Pk_ext); 

/*       for (j = 0; j < nk_emu; j++) */
/* 	{ */
/*       	  Pk_emu[j] *= pow(cosmo.h,3.); */
/* 	  Pk_emu[j] /= (2*M_PI*M_PI); */
/* 	  //	  fprintf(stderr, "%e %e\n", k2h[j], P2h[j]); */
/* 	} */


      for (j = 0; j < nk_ext; j++)
	{
      	  Pk_ext[j] *= pow(cosmo.h,3.);
	  Pk_ext[j] /= (2*M_PI*M_PI);
	  //	  fprintf(stderr, "%e %e\n", k2h[j], P2h[j]);
	}


      // need pade approximation // FIX THIS!!!
      //      pade(nk_emu, k_emu, Pk_emu, extra_pts, k_ext, Pk_ext); 

      // calculate chi 
      double chi = angular_diameter_distance(cosmo, cosmo.redshift);

      hankel_transform_gamma_t(hankel, cosmo, nk_ext, k_ext, Pk_ext, info.ntheta, info.theta, chi*c, info.pk_table[i]);

    }


  free(k_emu); free(Pk_emu);
  free(k_ext); free(Pk_ext); 
  return(0); 
}

int tabulate_pk_integral_DeltaSigma(struct survey_info survey, struct info_for_lens_dist_integrand info, struct HOD_params HOD, struct cosmo_params cosmo, struct info_for_hankel hankel)
{
  int i,j; 


  // set up arrays for emulator 

  int nk_emu = 661; 
  double *k_emu = malloc(nk_emu*sizeof(double)); 
  double *Pk_emu = malloc(nk_emu*sizeof(double)); 

  int extra_pts = 51; // no of points to extend in either direction 
  //  int nk_ext = nk_emu+2*extra_pts-2;
  int nk_ext = 976;
  double *k_ext  = malloc(nk_ext*sizeof(double));
  double *Pk_ext = malloc(nk_ext*sizeof(double));

  double params[5]; 
  convert_params(HOD, params); 

  for (i= 0; i < nk_emu; i++)
    {
      k_emu[i] = pow(10.,logk[i]); 
      k_emu[i]  /= cosmo.h;
      //      fprintf(stderr, "%lf\n", k_emu[i]); 
    }

  for (i = 0; i < nk_ext; i++)
    {
      k_ext[i] = pow(10.,k_m_ext[i]); 
      k_ext[i] /= cosmo.h;
    }

  double k_min = 1e-3;
  double k_max = 1e8;

  // These are determined by the snapshots used by the emulator: 
  double emu_scale_factors[6] = {0.9999,  0.9083, 0.8051, 0.6974,  0.6086, 0.4985};// last scalefactor is actually a = 1, but can't evaluate at chi = 0;  

  for (i = 0; i < info.nchi; i++)
    {
      // calculate the power spectrum from emulator

      cosmo.redshift = 1/emu_scale_factors[i] - 1; 

      emu(params, cosmo.redshift, Pk_ext); 

/*       for (j = 0; j < nk_emu; j++) */
/* 	{ */
/*       	  Pk_emu[j] *= pow(cosmo.h,3.); */
/* 	  Pk_emu[j] /= (2*M_PI*M_PI); */
/* 	  //	  fprintf(stderr, "%e %e\n", k2h[j], P2h[j]); */
/* 	} */


      for (j = 0; j < nk_ext; j++)
	{
      	  Pk_ext[j] *= pow(cosmo.h,3.);
	  Pk_ext[j] /= (2*M_PI*M_PI);
	  //	  fprintf(stderr, "%e %e\n", k2h[j], P2h[j]);
	}


      // need pade approximation // FIX THIS!!!
      //      pade(nk_emu, k_emu, Pk_emu, extra_pts, k_ext, Pk_ext); 

      // calculate chi 
      double chi = angular_diameter_distance(cosmo, cosmo.redshift);

      hankel_transform_DeltaSigma(hankel, cosmo, nk_ext, k_ext, Pk_ext, info.ntheta, info.theta, chi*c, info.pk_table[i]);

      /* char filename[256];  */
      /* sprintf(filename, "pk_check_%d.dat", i);  */
      /* FILE *fp_check = fopen(filename, "w");  */

      /* /\* for (j = 0; j < nk_ext; j++) *\/ */
      /* /\* 	fprintf(fp_check, "%lf %lf \n", k_ext[i], Pk_ext[i]);  *\/ */

      /* for (j = 0; j < nk_ext; j++) */
      /* 	fprintf(fp_check, "%lf %lf \n", k_ext[j], Pk_ext[j]);  */

      /* fclose(fp_check);  */
    }


  free(k_emu); free(Pk_emu);
  free(k_ext); free(Pk_ext); 
  return(0); 
}


int tabulate_fk_integral(int nbins_l, double *z_l, double *chi_l, double *n_l, int nbins_s, double *z_s, double *chi_s, double *n_s, int nchi, double *table)
{
  // Precompute array of all parts that only depend on chi
  // This is W_lens(chi)*f(chi)/a(chi)
  int i,j;


  double *f_K = malloc((nbins_l+1)*sizeof(double)); 
  double *scalefactor = malloc((nbins_l+1)*sizeof(double)); 

  // Precompute array of f_K(chi): 

  for (i = 0; i < nbins_l; i++)
    {
      /* chi[i] = log10(chi_min) + i*((log10(chi_max) - log10(chi_min))/((double)nchi-1)); */
      /* chi[i] = pow(10.,chi[i]); */
      /* fprintf(stderr, "%lf\n", chi[i]);  */
      f_K[i] = integrate_source_dist(chi_l[i], nbins_s, z_s, chi_s, n_s); 
      if (f_K[i] < 0)
	f_K[i] =0;
      //      fprintf(stderr, "%lf %5.3e\n", chi_l[i], f_K[i]);
    }

  /* // Normalise values?  */
  for (i = 0; i < nbins_l; i++)
    {
  /*     f_K[i] /= f_K[0]; */
      f_K[i] *= chi_l[i];
    }


  // Precompute array of scalefactors (as a function of radial distance)
  // spline values of chi as a function of z
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nbins_l);

  gsl_spline_init(spline, chi_l, z_l, nbins_l);


  for (i=0; i < nbins_l; i++)
    {
      double z = gsl_spline_eval(spline, chi_l[i], acc); 
      scalefactor[i] = 1./(1.+z); 
    }


  // Putting everything together: 

  for (i = 0; i < nbins_l; i++)
    {
      table[i] = c*n_l[i]*f_K[i]/scalefactor[i];  // table[i]/c/c is the same as g_tomo
      //      table[i] = n_l[i]; 
    } 

  /* FILE *fp_check = fopen("f_k_check_integral.dat", "w"); */
  /* for (i = 0; i < nbins_l; i++) */
  /*   fprintf(fp_check, "%lf %lf %5.3e %5.3e\n", chi_l[i], scalefactor[i], f_K[i], table[i]); */

  /* fclose(fp_check); */

  free(f_K); free(scalefactor); 
  gsl_spline_free(spline); gsl_interp_accel_free(acc); 
  return(0); 
}

double lens_dist_integrand(double chi, void *params)
{
  // n_l is the distribution of lenses as a function of chi
  // f_chi is the lensing weight (precomputed)
  // pk_grid_int is a table of precomputed integrals over J_2 and P(k)
  double f;
  struct info_for_lens_dist_integrand info = *(struct info_for_lens_dist_integrand*) (params);

  int i = 0;


  // evaluate splines

  double h_interp;
  double pk_grid_interp; 

  if (chi < info.chi[0])
      f = 0;
  /* else if (chi > info.chi[info.npts-1]) */
  /*     f = 0; */

  /* if (chi < info.chi_table[0]) */
  /*     f = 0; */
  else if (chi > info.chi_table[info.nchi-1])
      f = 0;
  else
    {
      pk_grid_interp = gsl_spline_eval(info.spline_pkJ, chi, info.acc_pkJ); 
      h_interp = gsl_spline_eval(info.spline_f_K, chi, info.acc_f_K); 
      //      f = pk_grid_interp;

      if (h_interp > 0)
      	f = h_interp * pk_grid_interp;
      else
      	f = 0;
    }

  //  fprintf(stderr, "%5.3e %5.3e %5.3e\n", h_interp, pk_grid_interp, f); 


  return(f);
}

int z_to_comoving(int nbins, double *z, double *chi, struct cosmo_params cosmo)
{
  //  double c = 2997.92458/cosmo.h; // Mpc

  // convert n(z) to n(chi)
  int i; 

  double error;

  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (1000);
  

  gsl_function F;
  F.function = &comoving_int;
  F.params = &cosmo.Om;

  for (i = 0; i < nbins; i++)
    {
      gsl_integration_qags (&F, 0, z[i], 0, 1e-7, 1000,  w, &(chi[i]), &error); 
      //      chi[i] *= c;
    }

  gsl_integration_workspace_free (w);


  return(0); 
}

double angular_diameter_distance(struct cosmo_params cosmo, double redshift)
{
  // this is the comoving angular diameter distance
  //  double c = 2997.92458/cosmo.h; // Mpc
  double da; 

  double comoving_dist; 
  double error;

  gsl_integration_cquad_workspace *w 
    = gsl_integration_cquad_workspace_alloc(1000);

  gsl_function F;
  F.function = &comoving_int;
  F.params = &cosmo.Om;

  gsl_integration_cquad(&F, 0, redshift, 0, 1e-3, w, &comoving_dist, &error, 0);

  gsl_integration_cquad_workspace_free(w);

  /* da = c*comoving_dist; */
   da = comoving_dist;

  return(da); 
}; 

double comoving_int(double z, void *params)
{
  double comoving; 
  double Omega_m = *(double *)params; 
  double Omega_lambda = 1-Omega_m; 

  //fprintf(stderr, "%lf %lf \n", Omega_m, Omega_lambda); 
  comoving = (1.+z)*(1.+z)*(1.+z)*Omega_m+Omega_lambda;
  comoving = 1./sqrt(comoving); 

  return(comoving); 
}

int setup_hankel_transform(struct info_for_hankel *hankel) 
{
  // set up information for hankel transform
  (*hankel).max_sum = 500; 
  (*hankel).zero     = malloc((*hankel).max_sum*sizeof(double)); 
  (*hankel).phi     = malloc((*hankel).max_sum*sizeof(double)); 
  (*hankel).phidash = malloc((*hankel).max_sum*sizeof(double)); 
  (*hankel).weights = malloc((*hankel).max_sum*sizeof(double)); 
  (*hankel).besselterm = malloc((*hankel).max_sum*sizeof(double)); 


  // Pre-compute values of functions not dependent on theta
  int n; 

  for (n = 1; n < (*hankel).max_sum; n++)
    {
      int nzero = n;
      double h = 1e-3;


      (*hankel).zero[n] = gsl_sf_bessel_zero_Jnu (2,nzero)/M_PI; // zero of Bessel/PI
      double hr = h*(*hankel).zero[n];

      (*hankel).phi[n] = hr *tanh(0.5*M_PI * sinh(hr));
      (*hankel).phidash[n] = tanh(0.5*M_PI*sinh(hr))+0.5*M_PI*hr*cosh(hr)/cosh(M_PI*0.5*sinh(hr))/cosh(M_PI*0.5*sinh(hr));
      (*hankel).weights[n] = gsl_sf_bessel_Ynu(2,M_PI*(*hankel).zero[n])/gsl_sf_bessel_Jnu(3,M_PI*(*hankel).zero[n]);
      (*hankel).besselterm[n] = gsl_sf_bessel_Jnu(2,M_PI*(*hankel).phi[n]/h);
    }

  return(0); 
}

int hankel_transform_gamma_t(struct info_for_hankel hankel, struct cosmo_params cosmo, int nk, double *k, double *Pk, int nr, double *theta, double chi, double *out)
{
  int n; 
  int i;
  double *Pkint = malloc(nk*sizeof(double)); 
  
  double sigma_smooth = 0.1;

  for (i =0; i < nk; i++)
    Pkint[i] = k[i]*Pk[i];//*exp(-k[i]*k[i]*sigma_smooth*sigma_smooth);

  gsl_interp_accel *acc 
    = gsl_interp_accel_alloc ();
  gsl_spline *spline 
    = gsl_spline_alloc (gsl_interp_cspline, nk);

  gsl_spline_init (spline, k, Pkint, nk);

  int m; 
  double da = chi;
  //double da = 1; 


  for (m = 0; m < nr; m++)
    {
      /* theta[m] = log10(r_min) + m*((log10(r_max) - log10(r_min))/(double)nr); // r = 80 is the limit for k = 1;  */
      /* theta[m] = pow(10.,theta[m]);  */

      double rtheta = theta[m]*da*arcminutes_to_radians; 
      double sum = 0; 
      /* for (n=0; n < nk; n++) */
      /* 	Pk[n] /= r;  */

      for (n = 1, sum=0; n < hankel.max_sum; n++)
	{
	  double nu = 2; 
	  int nzero = n; 
	  double zero = hankel.zero[n];
 	  double h = 1e-3; 
	  double hr = h*zero;
	  /* phi[n] = hr *tanh(0.5*M_PI * sinh(hr)); */
	  /* phidash[n] = tanh(0.5*M_PI*sinh(hr))+0.5*M_PI*hr*cosh(hr)/cosh(M_PI*0.5*sinh(hr))/cosh(M_PI*0.5*sinh(hr)); */

	  //	  double weights = gsl_sf_bessel_Ynu(2,M_PI*zero)/gsl_sf_bessel_Jnu(3,M_PI*zero);
	  /* double weights = -1.*sqrt(2./(M_PI*M_PI*zero))*cos(M_PI*zero)/(gsl_sf_bessel_Jnu(3,M_PI*zero)); */
	  double x = M_PI*hankel.phi[n]/h;
	  //double pkterm   = pow(x,1.5)*gsl_spline_eval (spline, x/r, acc);
	  double pkterm   = gsl_spline_eval (spline, x/rtheta, acc);
	  /* double besselterm = gsl_sf_bessel_Jnu(2,x); */
	  sum += M_PI*hankel.weights[n]*pkterm*hankel.besselterm[n]*hankel.phidash[n];
	  /* fprintf(stderr, "%d %f %f %f %f %f %f \n", n, sum, pkterm, besselterm, weights,  phi[n], derivs); */
	  /* fprintf(stderr, "%d %f %f %f %f \n", n, x, x/r, phi[n],  phidash[n] ); */
	}
      //      double constant = sqrt(0.5*M_PI)/(2.*M_PI*M_PI)/pow(r,1.5);
      out[m] = sum/rtheta; 


    }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free(Pkint); 


  return(0); 
}

int hankel_transform_DeltaSigma(struct info_for_hankel hankel, struct cosmo_params cosmo, int nk, double *k, double *Pk, int nr, double *theta, double chi, double *out)
{
  int n; 
  int i;
  double *Pkint = malloc(nk*sizeof(double)); 
  
  double sigma_smooth = 0.1;

  for (i =0; i < nk; i++)
    Pkint[i] = k[i]*Pk[i];//*exp(-k[i]*k[i]*sigma_smooth*sigma_smooth);

  gsl_interp_accel *acc 
    = gsl_interp_accel_alloc ();
  gsl_spline *spline 
    = gsl_spline_alloc (gsl_interp_cspline, nk);

  gsl_spline_init (spline, k, Pkint, nk);

  int m; 
  double da = chi;
  //double da = 1; 


  for (m = 0; m < nr; m++)
    {
      /* theta[m] = log10(r_min) + m*((log10(r_max) - log10(r_min))/(double)nr); // r = 80 is the limit for k = 1;  */
      /* theta[m] = pow(10.,theta[m]);  */
      double rtheta = theta[m];//theta[m]*da*arcminutes_to_radians; 
      double sum = 0; 
      /* for (n=0; n < nk; n++) */
      /* 	Pk[n] /= r;  */

      for (n = 1, sum=0; n < hankel.max_sum; n++)
	{
	  double nu = 2; 
	  int nzero = n; 
	  double zero = hankel.zero[n];
 	  double h = 1e-3; 
	  double hr = h*zero;
	  /* phi[n] = hr *tanh(0.5*M_PI * sinh(hr)); */
	  /* phidash[n] = tanh(0.5*M_PI*sinh(hr))+0.5*M_PI*hr*cosh(hr)/cosh(M_PI*0.5*sinh(hr))/cosh(M_PI*0.5*sinh(hr)); */

	  //	  double weights = gsl_sf_bessel_Ynu(2,M_PI*zero)/gsl_sf_bessel_Jnu(3,M_PI*zero);
	  /* double weights = -1.*sqrt(2./(M_PI*M_PI*zero))*cos(M_PI*zero)/(gsl_sf_bessel_Jnu(3,M_PI*zero)); */
	  double x = M_PI*hankel.phi[n]/h;
	  //double pkterm   = pow(x,1.5)*gsl_spline_eval (spline, x/r, acc);
	  double pkterm   = gsl_spline_eval (spline, x/rtheta, acc);
	  /* double besselterm = gsl_sf_bessel_Jnu(2,x); */
	  sum += M_PI*hankel.weights[n]*pkterm*hankel.besselterm[n]*hankel.phidash[n];
	  /* fprintf(stderr, "%d %f %f %f %f %f %f \n", n, sum, pkterm, besselterm, weights,  phi[n], derivs); */
	  /* fprintf(stderr, "%d %f %f %f %f \n", n, x, x/r, phi[n],  phidash[n] ); */
	}
      //      double constant = sqrt(0.5*M_PI)/(2.*M_PI*M_PI)/pow(r,1.5);
      out[m] = sum/rtheta; 


    }
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free(Pkint); 


  return(0); 
}



int convert_params(struct HOD_params HOD, double params[])
{
  // Emulator ranges
  float min_design[5] = {12.9, 13.5, 0.5, 0.5, 0.5};
  float max_design[5] = {14.0, 15.0, 1.2, 1.5, 1.5};


  params[0] = (HOD.Mcut  - min_design[0])/(max_design[0] - min_design[0]); 
  params[1] = (HOD.M1    - min_design[1])/(max_design[1] - min_design[1]); 
  params[2] = (HOD.sigma - min_design[2])/(max_design[2] - min_design[2]); 
  params[3] = (HOD.kappa - min_design[3])/(max_design[3] - min_design[3]); 
  params[4] = (HOD.alpha - min_design[4])/(max_design[4] - min_design[4]); 


  return(0); 
}
