struct info_for_hankel
{
  int max_sum; 
  double *zero;
  double *phi;
  double *phidash;
  double *weights;
  double *besselterm;
}; 

struct cosmo_params
{
  double Om; 
  double Omb;
  double w; 
  double ns; 
  double h; 
  double sigma8;
  double redshift;
  double growth_factor; 
}; 

struct HOD_params
{
  float Mcut;
  float M1;
  float sigma; 
  float alpha;
  float kappa;
}; 




struct survey_info
{
  int nbins_s; 
  double *z_s; 
  double *chi_s; 
  double *n_z_s; 
  double z_mean_s; 
  int nbins_l; 
  double *z_l; 
  double *chi_l; 
  double *n_z_l; 
  double z_mean_l; 
}; 

struct powerspectra
{
  int nk; 
  double *k; 
  double *Pk; 
}; 

struct info_for_source_dist_integrand
{
  int npts; 
  double *integrand; 
  double chi; 
  double z_l; 
}; 

struct info_for_lens_dist_integrand
{
  double *theta; // tracks which value of theta we are up to
  int npts; // number of bins for the lens distribution
  int nchi; // total number of points evaluated for the radial distance
  int ntheta; // total number of points evaluated for the angular scale
  double last; // last value of chi at which source distribution is non-zero
  double z_l;
  double *chi;
  double *chi_table; // this is for the power spectrum only
  double **pk_table; // holds results of int dk k P(k) J2
  double *f_chi_table; // hols results of collective factor n_l*f_K/a
  double *n_l; // lens distribution
  gsl_spline *spline_f_K; // spline for n_l*f_K/a
  gsl_interp_accel *acc_f_K; 
  gsl_spline *spline_pkJ; // spline for int dk k Pk J2
  gsl_interp_accel *acc_pkJ; 

}; 

struct info_for_sigma_crit_integrand
{
  struct cosmo_params cosmo; 
  struct survey_info survey; 
}; 
