/* Fitting Formulae for CDM + Baryon + Massive Neutrino (MDM) cosmologies. */
/* Daniel J. Eisenstein & Wayne Hu, Institute for Advanced Study */

/* There are two primary routines here, one to set the cosmology, the
other to construct the transfer function for a single wavenumber k. 
You should call the former once (per cosmology) and the latter as 
many times as you want. */

/* TFmdm_set_cosm() -- User passes all the cosmological parameters as
	arguments; the routine sets up all of the scalar quantites needed 
	computation of the fitting formula.  The input parameters are: 
	1) omega_matter -- Density of CDM, baryons, and massive neutrinos,
				in units of the critical density. 
	2) omega_baryon -- Density of baryons, in units of critical. 
	3) omega_hdm    -- Density of massive neutrinos, in units of critical 
	4) degen_hdm    -- (Int) Number of degenerate massive neutrino species 
	5) omega_lambda -- Cosmological constant 
	6) hubble       -- Hubble constant, in units of 100 km/s/Mpc 
	7) redshift     -- The redshift at which to evaluate */

/* TFmdm_onek_mpc() -- User passes a single wavenumber, in units of Mpc^-1.
	Routine returns the transfer function from the Eisenstein & Hu
	fitting formula, based on the cosmology currently held in the 
	internal variables.  The routine returns T_cb (the CDM+Baryon
	density-weighted transfer function), although T_cbn (the CDM+
	Baryon+Neutrino density-weighted transfer function) is stored
	in the global variable tf_cbnu. */

/* We also supply TFmdm_onek_hmpc(), which is identical to the previous
	routine, but takes the wavenumber in units of h Mpc^-1. */

/* We hold the internal scalar quantities in global variables, so that
the user may access them in an external program, via "extern" declarations. */

/* Please note that all internal length scales are in Mpc, not h^-1 Mpc! */

/* -------------------------- Prototypes ----------------------------- */

//int TFmdm_set_cosm(double omega_matter, double omega_baryon, double omega_hdm,
//	double degen_hdm, double omega_lambda, double hubble, double redshift);
//double TFmdm_onek_mpc(double kk);
//double TFmdm_onek_hmpc(double kk);

#include <math.h>
#include <stdio.h>
#include <cosmo.h>

/* Convenience from Numerical Recipes in C, 2nd edition */
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

/* ------------------------- Global Variables ------------------------ */

/* The following are set in TFmdm_onek_mpc() */
static double
	gamma_eff,	/* Effective \Gamma */
	growth_cb,	/* Growth factor for CDM+Baryon perturbations */
	growth_cbnu,	/* Growth factor for CDM+Baryon+Neutrino pert. */
	max_fs_correction,  /* Correction near maximal free streaming */
	qq,		/* Wavenumber rescaled by \Gamma */
	qq_eff,		/* Wavenumber rescaled by effective Gamma */
	qq_nu,		/* Wavenumber compared to maximal free streaming */
	tf_master,	/* Master TF */
	tf_sup,		/* Suppressed TF */
	y_freestream; 	/* The epoch of free-streaming for a given scale */

/* Finally, TFmdm_onek_mpc() and TFmdm_onek_hmpc() give their answers as */
static double   tf_cb,		/* The transfer function for density-weighted
			CDM + Baryon perturbations. */
	tf_cbnu;	/* The transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */

/* By default, these functions return tf_cb */

// This sets the redshift dependent parameters so that this can be done separately
//short TFmdm_set_cosm_change_z(double omega_matter, double omega_lambda, double redshift){
short COSMOLOGY::TFmdm_set_cosm_change_z(double redshift){
	double omega_denom;
	short qwarn = 0;
	struct EHparameters *ehparam = ehparam;

	if (redshift<=-1.0) {
		fprintf(stderr,"TFmdm_set_cosm(): Redshift < -1 is illegal.\n");
		/*exit(1);*/
	} else if (redshift>99.0) {
		fprintf(stderr,
		  "TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n");
		qwarn = 1;
	}

	// the below are dependent on redshift
	omega_denom = Oml+SQR(1.0+redshift)*(omega_curv + Omo*(1.0+redshift));
	omega_lambda_z = Oml/omega_denom;
	omega_matter_z = Omo*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
	growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/
			(pow(omega_matter_z,4.0/7.0)-omega_lambda_z+
		    (1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
	growth_to_z0 = z_equality*2.5*Omo/(pow(Omo,4.0/7.0)
		    -Oml + (1.0+Omo/2.0)*(1.0+Oml/70.0));
	growth_to_z0 = growth_k0/growth_to_z0;

	return qwarn;
}
/* ------------------------- TFmdm_set_cosm() ------------------------ */
//int TFmdm_set_cosm(double omega_matter, double omega_baryon, double omega_hdm,
//	double degen_hdm, double omega_lambda, double hubble)
short COSMOLOGY::TFmdm_set_cosm()
/* This routine takes cosmological parameters and a redshift and sets up
all the internal scalar quantities needed to compute the transfer function. */
/* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
				in units of the critical density. */
/* 	  omega_baryon -- Density of baryons, in units of critical. */
/* 	  omega_hdm    -- Density of massive neutrinos, in units of critical */
/* 	  degen_hdm    -- (Int) Number of degenerate massive neutrino species */
/*    omega_lambda -- Cosmological constant */
/* 	  hubble       -- Hubble constant, in units of 100 km/s/Mpc */
/*        redshift     -- The redshift at which to evaluate */
/* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
	sets many global variables for use in TFmdm_onek_mpc() */
{
    double z_drag_b1, z_drag_b2;
    int qwarn = 0;

  //  printf("omega_matter=%f omega_baryon=%f omega_hdm=%f \n"
  //  		,omega_matter,omega_baryon,omega_hdm);
  // 	printf("degen_hdm=%f omega_lambda=%e hubble=%e redshift=%e\n",
  //  		degen_hdm,omega_lambda,hubble,redshift);

    theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */

    /* Look for strange input */
    if (Omb<0.0) {
	fprintf(stderr,
	  "TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n");
	qwarn = 1;
    }
    if (Omnu<0.0) {
	fprintf(stderr,
	  "TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n");
	qwarn = 1;
    }
    if (h<=0.0) {
	fprintf(stderr,"TFmdm_set_cosm(): Negative Hubble constant illegal.\n");
	//exit(1);  /* Can't recover */
    } else if (h>2.0) {
	fprintf(stderr,"TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n");
	qwarn = 1;
    }

    if (Nnu<=0.0) Nnu=1;
    num_degen_hdm = (double) Nnu;
 	/* Have to save this for TFmdm_onek_mpc() */
    /* This routine would crash if baryons or neutrinos were zero, 
	so don't allow that */
    if (Omb<=0) Omb=1e-5;
    if (Omnu<=0) Omnu=1e-5;

    omega_curv = 1.0-Omo-Oml;
    omhh = Omo*SQR(h);
    obhh = Omb*SQR(h);
    onhh = Omnu*SQR(h);
    f_baryon = Omb/Omo;
    f_hdm = Omnu/Omo;
    f_cdm = 1.0-f_baryon-f_hdm;
    f_cb = f_cdm+f_baryon;
    f_bnu = f_baryon+f_hdm;

    /* Compute the equality scale. */
    z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));	/* Actually 1+z_eq */
    k_equality = 0.0746*omhh/SQR(theta_cmb);

    /* Compute the drag epoch and sound horizon */
    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*
		(1.0+z_drag_b1*pow(obhh,z_drag_b2));
    y_drag = z_equality/(1.0+z_drag);

    sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));

    /* Set up for the free-streaming & infall growth function */
    p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
    p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

    /* Compute small-scale suppression */
    alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*
	pow(1+y_drag,p_cb-p_c)*
	(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/
	(1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))*
	(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    alpha_gamma = sqrt(alpha_nu);
    beta_c = 1/(1-0.949*f_bnu);
    /* Done setting scalar variables */

   return qwarn;
}

/* ---------------------------- TFmdm_onek_mpc() ---------------------- */

double COSMOLOGY::TFmdm_onek_mpc(double kk)
/* Given a wavenumber in Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
/* Input: kk -- Wavenumber in Mpc^-1 */
/* Output: The following are set as global variables:
	growth_cb -- the transfer function for density-weighted
			CDM + Baryon perturbations. 
 	growth_cbnu -- the transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */
/* The function returns growth_cb */
{
    double tf_sup_L, tf_sup_C;
    double temp1, temp2;

    qq = kk/omhh*SQR(theta_cmb);

    /* Compute the scale-dependent growth functions */
    y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*
		SQR(num_degen_hdm*qq/f_hdm);
    temp1 = pow(growth_k0, 1.0-p_cb);
    temp2 = pow(growth_k0/(1+y_freestream),0.7);
    growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
    growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;

    /* Compute the master function */
    gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/
		(1+SQR(SQR(kk*sound_horizon_fit*0.43))));
//    printf("gamma_eff=%e\n",gamma_eff/omhh);
//    exit(0);

    qq_eff = qq*omhh/gamma_eff;

    tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
    tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
    tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));

    /******** modifications made by Kiakotou et al 2008 ***************/
    /*qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);*/
    qq_nu = 4.9*qq*sqrt(num_degen_hdm/f_hdm);

    /*max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/
      (pow(qq_nu,-1.6)+pow(qq_nu,0.8)); */
    max_fs_correction = 1+1.2*pow(f_hdm,0.68)*pow(num_degen_hdm,0.3+0.6*f_hdm)/
		(pow(qq_nu,-1.3)+pow(qq_nu,0.45));

    /******************************************************************/
    tf_master = tf_sup*max_fs_correction;

    /* Now compute the CDM+HDM+baryon transfer functions */
    tf_cb = tf_master*growth_cb/growth_k0;
    tf_cbnu = tf_master*growth_cbnu/growth_k0;
    return tf_cbnu;
}

/* ---------------------------- TFmdm_onek_hmpc() ---------------------- */

double COSMOLOGY::TFmdm_onek_hmpc(double kk)
/* Given a wavenumber in h Mpc^-1, return the transfer function for the
cosmology held in the global variables. */
/* Input: kk -- Wavenumber in h Mpc^-1 */
/* Output: The following are set as global variables:
	growth_cb -- the transfer function for density-weighted
			CDM + Baryon perturbations. 
 	growth_cbnu -- the transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */
/* The function returns growth_cb */
{
    return TFmdm_onek_mpc(kk*h);
}

#undef SQR(a)
