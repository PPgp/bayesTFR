#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*double test_rnorm(double *mu, double *sigma, double *t){
	double m, s;
	m = *mu;
	s = *sigma;

	GetRNGstate();
	t[0] = rnorm(m, s);
	PutRNGstate();
	GetRNGstate();

	PutRNGstate();
}*/

double sum(double *x, int dim) {
	double s;
	int i;
	s = 0.0;
	for (i=0; i<dim; ++i) s+=x[i];
	return(s);
}


void log_cond_Triangle_c4_trans(double *Triangle_c4_trans, double *eps_T, double *sd_eps_T, 
								double *mean_eps_T, int *dim_eps_T, double *Triangle4, double *delta4, 
								double *log_cond){
	int dimepsT, i;
	double s;
	
	dimepsT = *dim_eps_T;
	s = 0.0;
	for (i=0; i<dimepsT; i++) {
		s = s + dnorm(eps_T[i], mean_eps_T[i], sd_eps_T[i], 1);
	}
	*log_cond = -1/(2*pow(*delta4, 2.0)) * pow(*Triangle_c4_trans-(*Triangle4), 2.0) + s;
	/*Rprintf("\nlog_cond=%f", *log_cond);*/
}

void doDLcurve(double *DLpar, double *tfr, double *p1, double *p2, int *dim_tfr, 
				double *dl_values){
	double t_mid1, t_mid3, dl55, tmp1, tmp2;
	int i;
	
    t_mid1 = DLpar[3] + DLpar[2] + DLpar[1] + 0.5 * DLpar[0];
    t_mid3 = DLpar[3] + 0.5 * DLpar[2];
    dl55 = 5 * DLpar[4];    
    tmp1 = -log(pow(*p1, 2.0))/DLpar[0];
    tmp2 = -log(pow(*p2, 2.0))/DLpar[2];
    for (i=0; i< (*dim_tfr); i++){
    	if(tfr[i] <= 1) {dl_values[i] = 0;}
    	else {
    		dl_values[i] = dl55 * (-1/(1 + exp(tmp1 * (tfr[i] - t_mid1))) + 
    					1/(1 + exp(tmp2 * (tfr[i] - t_mid3))));
    		if(dl_values[i] < 0) dl_values[i] = 0;
    	}
    }
}

void dnormtrunc(double *x, double *mu, double *sigma, 
		double low, double high, int dim_out, double *out){
	int i;
	for (i=0; i< dim_out; i++) {
		if(x[i] < low || x[i] > high) out[i] = 0;
		else 
  		out[i] = dnorm(x[i],mu[i],sigma[i], 0)/(pnorm(high,mu[i],sigma[i],1,0)-pnorm(low,mu[i],sigma[i],1,0));
  	}
	return;
}



void dologdensityAR1(double *x, double *mu, double *sigma, 
			double *sigma_eps, double *rho_c, double *tfr, int *ltfr, double *low, double *high, double *logdens) {
	double dens[*ltfr-1], dnt[1];
	double s;
	int i;

	s = 0;
	for (i=0; i< (*ltfr - 1); i++){
		dens[i] = dnorm(tfr[i+1], *x + (*rho_c * (tfr[i]- *x)), *sigma_eps, 0);
		/*Rprintf("\n%f, %f %f %f", dens[i], dct[i], dl[i], loess_sd[i]);*/
		if(dens[i] < 1e-100) dens[i] = 1e-100;
		s = s+ log(dens[i]);
	}
	dnormtrunc(x, mu, sigma, *low, *high, 1, dnt);
	logdens[0] = s + log(dnt[0]);
	/*Rprintf("\ns=%f dnt=%f res = %f", s, dnt[0], logdens[0]);*/
	return;
}

