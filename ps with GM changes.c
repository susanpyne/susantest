#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/types.h>
#include "ps.h"
#include "lib.h"

extern cospar 		*cosmopar;
extern zdistr_par	*bg_zdistrpar;

double s1glob,s2glob;

extern double thetaglob;
extern double theta_min;
extern double theta_max;
extern double xi_plus_zero;

double thetaglob;
double theta_min=4.8e-6;//3e-5;
double theta_max=0.25;   //0.1;
double xi_plus_zero;
/*  global variable thetaglob is used here */
 

/*************************************************************************************/
/*  Functions for shear correlation functions                                                                                */
/*************************************************************************************/
  
 
double int_for_xi_plus(double x, void *args)
{
   //printf("%g\n",x*P_2(x/thetaglob)*ps_bessj0(x));
	return x*P_2(x/thetaglob)*ps_bessj0(x);
}

 
/*  global variable thetaglob is used here */
double int_for_xi_minus(double x, void *args)
{
	return x*P_2(x/thetaglob)*ps_bessj(4, x);
}


/* for theta=0 */
double int_for_xi_0(double s, void *args)
{
	return s*P_2(s);
}
  

/* integral for xi_plus */
double int_over_p2_j0(double theta)
{
 	double f, g, res, lower, upper;
	double RMAX, step;

	RMAX = 20.*theta*s_max;
	if (RMAX <= twopi)
		 RMAX = twopi+epsilon2; 
	if (theta < epsilon2)
	{
		f	=	int_GSL_integrate_qag(int_for_xi_0,NULL,0.0,1e8,NULL,2048)/twopi;
   	//f = ps_qromb(int_for_xi_0, 0.0, 1.e8)/twopi;
      return f;
	}
   
	thetaglob = theta;
	step = pi;
	lower = 2.35723;
	upper = step;
	f	=	int_GSL_integrate_qag(int_for_xi_plus,NULL,0.0,lower,NULL,2048);
	//f = ps_qromb1(int_for_xi_plus, 0.0, lower);
   
	while (1)
	{
		g = int_GSL_integrate_qag(int_for_xi_plus,NULL,lower,upper,NULL,2048);
		//g=ps_qromb1(int_for_xi_plus, lower, upper);
		if (fabs(g/f) < 1.e-6) 
			break;

		f += g;
		lower = upper;
		upper += step;
	}
	res = f/(theta*theta)/twopi;
	return res;
}


/* integral for xi_minus */
double int_over_p2_j4(double theta)
{
	double RMAX = 20.*theta*1.e6;
	double lower, upper, f, step, g, res;
	thetaglob = theta;

	if (theta<0.004) RMAX = 5.*0.004*1.e7;
	if (RMAX<twopi) { RMAX = twopi+epsilon2; }
	step = pi;
	lower = 2.35723;
	upper = step;
	f	=	int_GSL_integrate_qag(int_for_xi_minus,NULL,0.0,lower,NULL,2048);
	//f     = ps_qromb1(int_for_xi_minus, 0.0, lower);
	while (1)
	{
		g = int_GSL_integrate_qag(int_for_xi_minus,NULL,lower,upper,NULL,2048);
		//g = ps_qromb1(int_for_xi_minus, lower, upper);
		if (fabs(g/f) < 1.e-6) {
			/* printf("- %e %e\n", theta, upper); */
			break;
		}
		f += g;
		lower = upper;
		upper += step;
	}
	res = f/(twopi*(theta*theta));
	return res;
}


double xi(int pm, double theta)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double NONLINEAR = -42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table_p[N_theta];
	static double table_m[N_theta];
	static double dtheta = .0, logthetamin = 0., logthetamax = 0.;
	double res, thetalog, t;
	int i;
	FILE *F;
	char *xi_name;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || NONLINEAR != cosmopar->nonlinear || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
		logthetamin = log(theta_min);
		logthetamax = log(theta_max);
		dtheta = (logthetamax - logthetamin)/(N_theta - 1.);
      //printf("Func Xi: dtheta=%f\n",dtheta);
		xi_name = (char*)malloc(200*sizeof(char));
		sprintf(xi_name, "th/xi/xi-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f", cosmopar->omm, cosmopar->omv, cosmopar->Gamma, cosmopar->sigma8, cosmopar->n, bg_zdistrpar->beta, bg_zdistrpar->z0);
		/* read xi from file */
		if ((F=fopen(xi_name, "r")))
		{
			fscanf(F, "0.0 %lf\n", &xi_plus_zero);
			thetalog = logthetamin;
			i = 0;
			while (fscanf(F, "%lf %lf %lf\n", &t, table_p+i, table_m+i) != EOF)
			{
				if (fabs(t - exp(thetalog)) > epsilon3)
				{
					printf("xi: t:%g, tlog: %g\n",t,exp(thetalog));
					printf("Inconsistency with xi file %s!\n",xi_name); exit(1);
				}
				thetalog += dtheta;
				i++;
			}
			fileclose(F);
		}
		else
		{
			/* calculate xi and write to disk */
			F = fileopen(xi_name, "w");
			thetalog = logthetamin;
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
				t = exp(thetalog);
				table_p[i] = int_over_p2_j0(t);
				table_m[i] = int_over_p2_j4(t);
				printf("Th. cf: %1.1f%% done.\r",1.0*i/(1.0*N_theta)*100.); fflush(stdout);
			}
			xi_plus_zero = int_over_p2_j0(0.0);
			printf("Th. cf: 100.0%% done.\n");
			fprintf(F, "0.0 %5.20e\n", xi_plus_zero);
			thetalog = logthetamin;
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
				fprintf(F, "%5.10e %5.20e %5.20e\n", exp(thetalog), table_p[i], table_m[i]);
			}
			fileclose(F);
			free((void*)xi_name);
		}
	}
	/* interpolate xi */
	if (theta>=theta_min)
	{
		if (theta>theta_max)
		{
			printf("theta=%e larger than theta_max=%e in xi\n", theta, theta_max);
			error("theta too large in xi");
		}
		thetalog = log(theta);
		t = (thetalog-logthetamin)/dtheta;
		i = (int)(floor(t));
		if (pm == +1)
		{
			res = (t-i)*(table_p[i+1] - table_p[i]) + table_p[i];
		}
		else
		{
			res = (t-i)*(table_m[i+1] - table_m[i]) + table_m[i];
		}
		return res;
	}
	else if (theta < epsilon2)
	{
		if (pm == +1) { return xi_plus_zero; }
		else { return 0.0; }
	}
	else
	{
		t = theta/theta_min;
		if (pm == +1) { return t*(table_p[0]-xi_plus_zero) + xi_plus_zero; }
		else { return t*(table_m[0]); }
	}
}


double xi2(int pm, double theta,char *path)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double NONLINEAR = -42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table_p[N_theta];
	static double table_m[N_theta];
	static double dtheta = .0, logthetamin = 0., logthetamax = 0.;
	double res, thetalog, t;
	int i;
	FILE *F;
	char *xi_name;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || NONLINEAR != cosmopar->nonlinear || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
		logthetamin = log(theta_min);
		logthetamax = log(theta_max);
		dtheta = (logthetamax - logthetamin)/(N_theta - 1.);
      //printf("Func Xi: dtheta=%f\n",dtheta);
		xi_name = (char*)malloc(200*sizeof(char));
		sprintf(xi_name, "%s/xi-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f",path,cosmopar->omm, cosmopar->omv, cosmopar->Gamma, cosmopar->sigma8, cosmopar->n, bg_zdistrpar->beta, bg_zdistrpar->z0);
		/* read xi from file */
		if ((F=fopen(xi_name, "r")))
		{
			fscanf(F, "0.0 %lf\n", &xi_plus_zero);
			thetalog = logthetamin;
			i = 0;
			while (fscanf(F, "%lf %lf %lf\n", &t, table_p+i, table_m+i) != EOF)
			{
				if (fabs(t - exp(thetalog)) > epsilon3)
				{
					printf("xi: t:%g, tlog: %g\n",t,exp(thetalog));
					printf("Inconsistency with xi file %s!\n",xi_name); exit(1);
				}
				thetalog += dtheta;
				i++;
			}
			fileclose(F);
		}
		else
		{
			/* calculate xi and write to disk */
			F = fileopen(xi_name, "w");
			thetalog = logthetamin;
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
				t = exp(thetalog);
				table_p[i] = int_over_p2_j0(t);
				table_m[i] = int_over_p2_j4(t);
				printf("Th. cf: %1.1f%% done.\r",1.0*i/(1.0*N_theta)*100.); fflush(stdout);
			}
			xi_plus_zero = int_over_p2_j0(0.0);
			printf("Th. cf: 100.0%% done.\n");
			fprintf(F, "0.0 %5.20e\n", xi_plus_zero);
			thetalog = logthetamin;
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
				fprintf(F, "%5.10e %5.20e %5.20e\n", exp(thetalog), table_p[i], table_m[i]);
			}
			fileclose(F);
			free((void*)xi_name);
		}
	}
	/* interpolate xi */
	if (theta>=theta_min)
	{
		if (theta>theta_max)
		{
			printf("theta=%e larger than theta_max=%e in xi\n", theta, theta_max);
			error("theta too large in xi");
		}
		thetalog = log(theta);
		t = (thetalog-logthetamin)/dtheta;
		i = (int)(floor(t));
		if (pm == +1)
		{
			res = (t-i)*(table_p[i+1] - table_p[i]) + table_p[i];
		}
		else
		{
			res = (t-i)*(table_m[i+1] - table_m[i]) + table_m[i];
		}
		return res;
	}
	else if (theta < epsilon2)
	{
		if (pm == +1) { return xi_plus_zero; }
		else { return 0.0; }
	}
	else
	{
		t = theta/theta_min;
		if (pm == +1) { return t*(table_p[0]-xi_plus_zero) + xi_plus_zero; }
		else { return t*(table_m[0]); }
	}
}



/*************************************************************************************/
/*  Tomography: Functions for shear cross-correlation functions                      */
/*************************************************************************************/

 
double int_for_xi_plusx(double x, void *args)
{
	int b1,b2;
	double *arg	=	(double*)args;
	b1=(int)arg[0];
	b2=(int)arg[1];

	return x*P_2x(x/thetaglob,b1,b2)*ps_bessj0(x);
}
 

double int_for_xi_minusx(double x, void *args)
{
	int b1,b2;
	double *arg	=	(double*)args;
	b1=(int)arg[0];
	b2=(int)arg[1];
  
	return x*P_2x(x/thetaglob,b1,b2)*ps_bessj(4, x);
}


/* for theta=0 */
double int_for_xi_0x(double s, void *args)
{
	double *arg	=	(double*)args;
	
	int b1,b2;
	b1=(int)arg[0];
	b2=(int)arg[1];
  
	return s*P_2x(s,b1,b2);
}
  

/* integral for xi_plus */
double int_over_p2_j0x(double theta, void *args)
{
   //printf("int_over_p2_j0\n");
	double f, g, res, lower, upper;
	double RMAX, step;
	double *arg	=	(double*)args;
   
	RMAX = 20.*theta*s_max;
	if (RMAX <= twopi) { RMAX = twopi+epsilon2; }
	if (theta < epsilon2)
	{
		f	=	int_GSL_integrate_qag(int_for_xi_0x,NULL,0.0,1e8,NULL,2048)/twopi;
		//f = int_qrombn(int_for_xi_0x, 0.0, 1.e8,arg,argc)/twopi;

		return f;
	}
   
	thetaglob = theta;
	step = pi;
	lower = 2.35723;
	upper = step;
	f = int_GSL_integrate_qag(int_for_xi_plusx,NULL,0.0,lower,NULL,2048);//int_qrombn(int_for_xi_plusx, 0.0, lower,arg,argc);

	while (1)
	{
      g = int_GSL_integrate_qag(int_for_xi_plusx,NULL,lower,upper,NULL,2048);//int_qrombn(int_for_xi_plusx, lower, upper,arg,argc);
		if (fabs(g/f) < 1.e-6) 
			break;
	
		f += g;
		lower = upper;
		upper += step;
	}
	res = f/(theta*theta)/twopi;
	return res;
}


/* integral for xi_minus */
double int_over_p2_j4x(double theta,void *args)
{
	double RMAX = 20.*theta*1.e6;
	double lower, upper, f, step, g, res;
	double *arg	=	(double*)args;
	thetaglob = theta;

	if (theta<0.004) RMAX = 5.*0.004*1.e7;
	if (RMAX<twopi) { RMAX = twopi+epsilon2; }
	step = pi;
	lower = 2.35723;
	upper = step;
	f     = int_GSL_integrate_qag(int_for_xi_minusx,NULL,0.0,lower,NULL,2048);//int_qrombn(int_for_xi_minusx, 0.0, lower,arg,argc);
	while (1) {
		g = int_GSL_integrate_qag(int_for_xi_minusx,NULL,lower,upper,NULL,2048);//int_qrombn(int_for_xi_minusx, lower, upper,arg,argc);
		if (fabs(g/f) < 1.e-6)
			break;

		f += g;
		lower = upper;
		upper += step;
	}
	res = f/(twopi*(theta*theta));
	return res;
}


double xix(int pm, double theta, int bin1, int bin2)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double NONLINEAR = -42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static double table_p[N_theta];
	static double table_m[N_theta];
	static double dtheta = .0, logthetamin = 0., logthetamax = 0.;
	double arg[2];
	double res, thetalog, t;
	int i;
	FILE *F;
	char *xi_name;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || NONLINEAR != cosmopar->nonlinear || NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		NBIN1 = bin1;
		NBIN2 = bin2;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
		logthetamin = log(theta_min);
		logthetamax = log(theta_max);
		dtheta = (logthetamax - logthetamin)/(N_theta - 1.);
      
		arg[0]=bin1;
		arg[1]=bin2;

      //printf("Func Xi: dtheta=%f\n",dtheta);
		xi_name = (char*)malloc(300*sizeof(char));
		sprintf(xi_name, "th/xi/tomo/xi_%d-%d_-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f", bin1, bin2, cosmopar->omm, cosmopar->omv, cosmopar->Gamma, cosmopar->sigma8, cosmopar->n, bg_zdistrpar->beta, bg_zdistrpar->z0);
		/* read xi from file */
		if ((F=fopen(xi_name, "r")))
		{
			fscanf(F, "0.0 %lf\n", &xi_plus_zero);
			thetalog = logthetamin;
			i = 0;
			while (fscanf(F, "%lf %lf %lf\n", &t, table_p+i, table_m+i) != EOF)
			{
				if (fabs(t - exp(thetalog)) > epsilon3)
				{
					printf("xi: t:%g, tlog: %g\n",t,exp(thetalog));
					printf("Inconsistency with xi file %s!\n",xi_name); exit(1);
				}
				thetalog += dtheta;
				i++;
			}
			fileclose(F);
		}
		else
		{
			/* calculate xi and write to disk */
			F = fileopen(xi_name, "w");
			thetalog = logthetamin;
	 //printf("172\n");
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
	    //printf("i=%d\n",i);
				t = exp(thetalog);
				table_p[i] = int_over_p2_j0x(t,arg);
				table_m[i] = int_over_p2_j4x(t,arg);
				printf("Th. cf: %1.1f%% done.\r",1.0*i/(1.0*N_theta)*100.); fflush(stdout);
	    //printf("theta: %g, p: %g, m: %g\n",t,table_p[i],table_m[i]);	
			}
			xi_plus_zero = int_over_p2_j0x(0.0,arg);
			printf("Th. cf: 100.0%% done.\n");
			fprintf(F, "0.0 %5.20e\n", xi_plus_zero);
			thetalog = logthetamin;
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
				fprintf(F, "%5.10e %5.20e %5.20e\n", exp(thetalog), table_p[i], table_m[i]);
			}
			fileclose(F);
			free((void*)xi_name);
		}
	}
	/* interpolate xi */
	if (theta>=theta_min)
	{
		if (theta>theta_max)
		{
			printf("theta=%e larger than theta_max=%e in xi\n", theta, theta_max);
			error("theta too large in xi");
		}
		thetalog = log(theta);
		t = (thetalog-logthetamin)/dtheta;
		i = (int)(floor(t));
		if (pm == +1)
		{
			res = (t-i)*(table_p[i+1] - table_p[i]) + table_p[i];
		}
		else
		{
			res = (t-i)*(table_m[i+1] - table_m[i]) + table_m[i];
		}
		return res;
	}
	else if (theta < epsilon2)
	{
		if (pm == +1) { return xi_plus_zero; }
		else { return 0.0; }
	}
	else
	{
		t = theta/theta_min;
		if (pm == +1) { return t*(table_p[0]-xi_plus_zero) + xi_plus_zero; }
		else { return t*(table_m[0]); }
	}
}


double xix2(int pm, double theta, int bin1, int bin2,char *path)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double NONLINEAR = -42;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table_p[N_theta];
	static double table_m[N_theta];
	static double dtheta = .0, logthetamin = 0., logthetamax = 0.;
	double arg[2];
	double res, thetalog, t;
	int i;
	FILE *F;
	char *xi_name;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || NONLINEAR != cosmopar->nonlinear || NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		NBIN1 = bin1;
		NBIN2 = bin2;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
		logthetamin = log(theta_min);
		logthetamax = log(theta_max);
		dtheta = (logthetamax - logthetamin)/(N_theta - 1.);

		arg[0]=bin1;
		arg[1]=bin2;
 
      //printf("Func Xi: dtheta=%f\n",dtheta);
		xi_name = (char*)malloc(300*sizeof(char));
		mkdir(path,0777);
		sprintf(xi_name, "%s/xi_%d-%d_-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f", path, bin1, bin2, cosmopar->omm, cosmopar->omv, cosmopar->Gamma, cosmopar->sigma8, cosmopar->n, bg_zdistrpar->beta, bg_zdistrpar->z0);
		/* read xi from file */
		if ((F=fopen(xi_name, "r")))
		{
			fscanf(F, "0.0 %lf\n", &xi_plus_zero);
			thetalog = logthetamin;
			i = 0;
			while (fscanf(F, "%lf %lf %lf\n", &t, table_p+i, table_m+i) != EOF)
			{
				if (fabs(t - exp(thetalog)) > epsilon3)
				{
					printf("xi: t:%g, tlog: %g\n",t,exp(thetalog));
					printf("Inconsistency with xi file %s!\n",xi_name); exit(1);
				}
				thetalog += dtheta;
				i++;
			}
			fileclose(F);
		}
		else
		{
			/* calculate xi and write to disk */
			F = fileopen(xi_name, "w");
			thetalog = logthetamin;
	 //printf("172\n");
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
	    //printf("i=%d\n",i);
				t = exp(thetalog);
				table_p[i] = int_over_p2_j0x(t,arg);
				table_m[i] = int_over_p2_j4x(t,arg);
				printf("Th. cf: %1.1f%% done.\r",1.0*i/(1.0*N_theta)*100.); fflush(stdout);
	    //printf("theta: %g, p: %g, m: %g\n",t,table_p[i],table_m[i]);	
			}
			xi_plus_zero = int_over_p2_j0x(0.0,arg);
			printf("Th. cf: 100.0%% done.\n");
			fprintf(F, "0.0 %5.20e\n", xi_plus_zero);
			thetalog = logthetamin;
			for (i=0; i<N_theta; i++, thetalog+=dtheta)
			{
				fprintf(F, "%5.10e %5.20e %5.20e\n", exp(thetalog), table_p[i], table_m[i]);
			}
			fileclose(F);
			free((void*)xi_name);
		}
	}
	/* interpolate xi */
	if (theta>=theta_min)
	{
		if (theta>theta_max)
		{
			printf("theta=%e larger than theta_max=%e in xi\n", theta, theta_max);
			error("theta too large in xi");
		}
		thetalog = log(theta);
		t = (thetalog-logthetamin)/dtheta;
		i = (int)(floor(t));
		if (pm == +1)
		{
			res = (t-i)*(table_p[i+1] - table_p[i]) + table_p[i];
		}
		else
		{
			res = (t-i)*(table_m[i+1] - table_m[i]) + table_m[i];
		}
		return res;
	}
	else if (theta < epsilon2)
	{
		if (pm == +1) { return xi_plus_zero; }
		else { return 0.0; }
	}
	else
	{
		t = theta/theta_min;
		if (pm == +1) { return t*(table_p[0]-xi_plus_zero) + xi_plus_zero; }
		else { return t*(table_m[0]); }
	}
}




/* ============================================================ *
 * Global variables needed in integrand functions int_for_p_2	*
 * and int_for_g						*
 * ============================================================ */
double sglob, aglob;


/* CPT 9 */
double da_dtau(double a)
{
	double res, det;
	if (a < epsilon3) error("Division by 0 in da_dtau");
	det = 1. + cosmopar->omm*(1./a - 1.) + OMV_MACRO(a)*(a*a - 1.);
	res = sqrt(det);
	return res;
}


double da_dtau_m3(double a, void *args)
{
	double res;
	if (a < epsilon3)
		return 0.;
	res = da_dtau(a);
	if (res < epsilon3)
		error("Division by 0 in da_dtau_m3");
	return 1./(res*res*res);
}


double D_plus(double a)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double da = 0.;
	static double table[N_a];
	double delta, aa, delta0;
	int  i;

	if (fabs(cosmopar->omm-OMEGA_M)>epsilon3 || fabs(cosmopar->omv-OMEGA_V)>epsilon3 || W0 != cosmopar->w0 || WA != cosmopar->wa) {

		delta0 = int_GSL_integrate_qag(da_dtau_m3,NULL,0.,1.,NULL,2048);;

		da = (1. - a_min)/(N_a-1.);
		aa = a_min;
		for (i = 0; i<N_a; i++, aa += da)
		{
			delta = int_GSL_integrate_qag(da_dtau_m3,NULL,0.,aa,NULL,2048);
			table[i] = 1./aa * da_dtau(aa) * delta / delta0;
		}
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
	}
	return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}




/* BBKS G3, see also PD2 51, 52. Returns k^n * T^2(k) */
double Tsqr(double k)  //Transferfunktn^2	// k in Hubble length-1
{
  double q, f1, f2;
  const double coverh=ckms/100.;
  q = k/(coverh*cosmopar->Gamma);	// k in Hubble length-1, k/3000 in hMpc-1
  f1 = log(1 + 2.34*q)/(2.34*q);
  f2 = 1 + q*(3.89 + q*(259.21 + q*(162.771336 + q*2027.16958081)));
  /* polynomial: 1 + 3.89q + (16.1q)^2 + (5.46q)^3 + (6.71q)^4 */
  //assert(finite(f1)&&finite(f2));
 
  return pow(k, cosmopar->n)*(f1*f1)/sqrt(f2);
}



/*Eisenstein &Hu no-wiggle version*/
/*
double Tsqr_EH(double k)		// from Martin White, approximation for s as in EH98, (26)
{
	double q, theta, ommh2, a, s, gamma, L0, C0;
	double tmp;
	double omegam, ombh2, hubble;

	// other input parameters 
	hubble = cosmopar->h;

	omegam = cosmopar->omm;
	ombh2 = cosmopar->omb* cosmopar->h *cosmopar->h ;

	if(cosmopar->omb == 0)
		ombh2 = 0.04 * cosmopar->h *cosmopar->h ;

	k *= 1./3000.;//(3.085678e24 / UnitLength_in_cm);	// convert to h/Mpc

	theta = 2.728 / 2.7;
	ommh2 = omegam * hubble * hubble;
	s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
	a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
			+ 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
	gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
	gamma *= omegam * hubble;
	q = k * theta * theta / gamma; 
	L0 = log(2. * exp(1.) + 1.8 * q);
	C0 = 14.2 + 731. / (1. + 62.5 * q);
	tmp = L0 / (L0 + C0 * q * q);
	return (tmp*tmp)*pow(k,cosmopar->n);
}
*/

double Tsqr_EH(double k)		// as in EH98
{
  double q, theta, thetasq, ommh2, a, s, gamma, L0, C0;
	double tmp;
	double omegam, ombh2, hubble;
  const double coverh=ckms/100.;

	hubble = cosmopar->h;

	omegam = cosmopar->omm;
	ombh2 = cosmopar->omb* cosmopar->h *cosmopar->h ;

	if(cosmopar->omb == 0)
		ombh2 = 0.04 * cosmopar->h *cosmopar->h ;

	k *= 1./coverh;//(3.085678e24 / UnitLength_in_cm);	// convert to h/Mpc    

	theta = 2.728 / 2.7;
	thetasq=theta*theta;
	ommh2 = omegam * hubble * hubble;

	// compute s
	double z_eq,z_d,R_eq,R_d,k_eq,b_1,b_2;

	k_eq=7.46e-2*ommh2/thetasq;
	b_1=0.313*pow(ommh2,-0.419)*(1.+0.607*pow(ommh2,0.674));
	b_2=0.238*pow(ommh2,0.223);
	z_eq=2.5e4*ommh2/(thetasq*thetasq);
	z_d=1291.*pow(ommh2,0.251)*(1.+b_1*pow(ombh2,b_2))/(1.+0.659*pow(ommh2,0.828));
	R_eq=3.15e4*ombh2/(z_eq*thetasq*thetasq);
	R_d=3.15e4*ombh2/(z_d*thetasq*thetasq);
	s = 2./(3.*k_eq)*sqrt(6./R_eq)*log((sqrt(1.+R_d)+sqrt(R_d+R_eq))/(1.+sqrt(R_eq)));

	a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
			+ 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
	gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k*hubble * s)));  //this k is 1/Mpc
	gamma *= omegam * hubble;
	q = k * thetasq / gamma; 
	L0 = log(2. * exp(1.) + 1.8 * q);
	C0 = 14.2 + 731. / (1. + 62.5 * q);
	tmp = L0 / (L0 + C0 * q * q);

	return (tmp*tmp)*pow(k,cosmopar->n);
}
 

// Eh transfer function fit with wiggles
double T_tilde(double k,double a, double b)
{
  double L0,C0,q,theta;

  theta = 2.728 / 2.7;
  q=k*theta*theta/(cosmopar->omm*cosmopar->h);
  L0=log(exp(1.)+1.8*b*q);
  C0=14.2/a+386./(1.+69.9*pow(q,1.08));
  return(L0/(L0+C0*q*q));
}

double Gfunc_EH98(double x)
{
  double sx=sqrt(1.+x);
  return(x*((2.+3.*x)*log((sx+1.)/(sx-1.))-6.*sx));
}


double Tsqr_EH_wiggle(double k)		// as in EH98, eqs. 11-15,17-24
{
  double theta, thetasq, ombh2, ommh2, fb, fc, T_cdm,T_baryon,T_total;
  const double coverh=ckms/100.;

  if(cosmopar->omb == 0) {
    printf("Error: set Omega_b>0 for EH transfer function with wiggles.\n");
    exit(-1);
  }

  ombh2 = cosmopar->omb* cosmopar->h *cosmopar->h;
  ommh2 = cosmopar->omm* cosmopar->h *cosmopar->h;

  k *= 1./coverh;   // convert to h/Mpc    

  theta = 2.728 / 2.7;
  thetasq=theta*theta;

  fb=cosmopar->omb/cosmopar->omm;
  fc=(cosmopar->omm-cosmopar->omb)/cosmopar->omm;

  // compute sound horizon s
  double z_eq,z_d,R_eq,R_d,k_eq,b_1,b_2,s;

  k_eq=7.46e-2*ommh2/thetasq/cosmopar->h;  // divide by h to get units [h/Mpc]
  b_1=0.313*pow(ommh2,-0.419)*(1.+0.607*pow(ommh2,0.674));
  b_2=0.238*pow(ommh2,0.223);
  z_eq=2.5e4*ommh2/(thetasq*thetasq);
  z_d=1291.*pow(ommh2,0.251)*(1.+b_1*pow(ombh2,b_2))/(1.+0.659*pow(ommh2,0.828));
  R_eq=3.15e4*ombh2/(z_eq*thetasq*thetasq);
  R_d=3.15e4*ombh2/(z_d*thetasq*thetasq);
  s = 2./(3.*k_eq)*sqrt(6./R_eq)*log((sqrt(1.+R_d)+sqrt(R_d+R_eq))/(1.+sqrt(R_eq)));  // in units because of k_eq [Mpc/h]

  // CDM TF parameters
  double a1,a2,alpha_c,b1,b2,beta_c,f;
  a1=pow(46.9*ommh2,0.670)*(1.+pow(32.1*ommh2,-0.532));
  a2=pow(12.0*ommh2,0.424)*(1.+pow(45.0*ommh2,-0.582));
  alpha_c=pow(a1,(-1.)*fb)*pow(a2,(-1.)*pow(fb,3.));
  b1=0.944/(1.+pow(458.*ommh2,-0.708));
  b2=pow(0.395*ommh2,-0.0266);
  beta_c=1./(1.+b1*(pow(fc,b2)-1.));
  f=1./(1.+pow(k*s/5.4,4.));
  T_cdm=f*T_tilde(k,1.,beta_c)+(1.-f)*T_tilde(k,alpha_c,beta_c);

  // baryon TF parameters
  double alpha_b,beta_b,beta_node,s_tilde,k_silk;
  alpha_b=2.07*k_eq*s*pow(1.+R_d,-0.75)*Gfunc_EH98((1.+z_eq)/(1.+z_d));
  beta_b=0.5+fb+(3.-2.*fb)*sqrt(pow(17.2*ommh2,2.)+1.);
  beta_node=8.41*pow(ommh2,0.435);
  s_tilde=s/pow(1.+pow(beta_node/(k*s),3.),1./3.);
  k_silk=1.6*pow(ombh2,0.52)*pow(ommh2,0.73)*(1.+pow(10.4*ommh2,-0.95))/cosmopar->h;  // divide by h to get units [h/Mpc]
  T_baryon=sin(k*s_tilde)/(k*s_tilde)*(T_tilde(k,1.,1.)/(1.+pow(k*s/5.2,2.))+alpha_b/(1.+pow(beta_b/(k*s),3.))*exp((-1.)*pow(k/k_silk,1.4)));
  
  T_total=fb*T_baryon+fc*T_cdm;
  return(T_total*T_total)*pow(k,cosmopar->n);
}
 

void	LineFit(double x1, double x2, double y1, double y2, double *m, double *b)
{
	*m	=	(y2-y1)/(x2-x1);
	*b	=	y1-*m*x1;

	return;
}


double Tsqr_tabulated(double k)
{
	static int		firststart=1,nbin=0;
	static double	*kk,*T, lm, lb, um, ub;
	double			kint,Tint;
	int				i,j;
	const double coverh=ckms/100.;

	k *= 1./coverh; /* convert to h/Mpc */
	
	if(firststart)
	{
		printf("  Using tabulateted tranfer function\n");
		
		kk	=	calloc(1,sizeof(double));
		T	=	calloc(1,sizeof(double));
		
		FILE *fp	=	fopen(cosmopar->tffile,"r");
		if(!fp)
		{
			printf("  Error: Couldn't open %s!\n",cosmopar->tffile);
			exit(1);
		}
		
		while(!feof(fp))
		{
			fscanf(fp,"%lg %lg\n",&(kint),&(Tint));

			if(kint>=cosmopar->tfmin && kint <= cosmopar->tfmax)
			{
				kk	=	realloc(kk,(nbin+1)*sizeof(double));
				T	=	realloc(T,(nbin+1)*sizeof(double));
				printf(" %g %g\n",kint,Tint);
				kk[nbin]	=	kint;
				T[nbin]	=	Tint;

				nbin++;
			}
		}
		
		if(nbin<10)
		{
			printf("  Error: less than 10 valid data points in transferfunc. file!\n");
			exit(1);
		}
		fclose(fp);

		printf("  %d bins read.\n",nbin);
		
		//Prepare extraploation
		double mint, bint;
		for(i=0;i<4;i++)
		{
			LineFit(log(kk[i]),log(kk[i+1]),log(T[i]),log(T[i+1]),&mint,&bint);
			lm	+=	mint/4.; lb	+=	bint/4.;

			LineFit(log(kk[nbin-i-2]),log(kk[nbin-i-3]),log(T[nbin-i-2]),log(T[nbin-i-3]),&mint,&bint);
			um	+=	mint/4.; ub	+=	bint/4.;
		}
		
		firststart = 0;
		printf("  Prep done: %g %g %g %g\n", lm,lb,um,ub);
	}

	double logk	=	log(k);
	
	//interpolate from table, if k in range
	if(k>cosmopar->tfmin && k<cosmopar->tfmax)
	{
		int		lower, upper;
		double	dlogk,lkl,lku,u;

		//find nearest smaller table entry by bisection
		lower=0; upper=nbin-1;
		do
		{
			if(k<=kk[(lower+upper)/2])
				upper	=	(lower+upper)/2;
			else
				lower = (lower+upper)/2;
		} while (upper-lower>1);

		//interpolate
		lkl	=	log(kk[lower]);
		lku	=	log(kk[lower+1]);
		dlogk	=	lku-lkl;
		u	=	(logk - lkl)/dlogk;

		return exp( (1-u)*log(T[lower]) + u*log(T[lower+1]) );
	}
	//else extrapolate
	else if(k<=cosmopar->tfmin)
	{
		return exp( lm*logk + lb );
	}
	else if(k>=cosmopar->tfmin)
	{
		return exp( um*logk + ub );
	}

	return 0;
}


double int_for_sigma_8(double k, void *args)
{
  double kR, res, x;
  const double coverh=ckms/100.;
	
  kR = k/(coverh/8.);	//previously: 375=(3000.h-1Mpc/Hubble radius)/(8h-1Mpc)
  x = (sin(kR) - kR*cos(kR))/(kR*kR*kR);//filter function

  if(cosmopar->transfer_tabulated) {
    res=k*k*Tsqr_tabulated(k)*x*x;
  }
  else {
    if (cosmopar->transfer_EH==1) res = k*k*Tsqr_EH(k)*x*x;  //EH
    else if (cosmopar->transfer_EH==2) res = k*k*Tsqr_EH_wiggle(k)*x*x;  //EHW
    else res = k*k*Tsqr(k)*x*x;  //BBKS/EBW
  }
  return res;
}


/* PD2 42, for a=1, so D+=1 */
double sigma_8_sqr()    	//for normalizing P_L
{
	static double N_SPEC = -42.;
	static double GAMMA  = -42.;
	static double res = -42.;
	double integral,keff;

	if (fabs(N_SPEC-cosmopar->n)>epsilon3 || fabs(GAMMA-cosmopar->Gamma)>epsilon3)
	{
          // full integral - method of choice!
	  integral = int_GSL_integrate_qag(int_for_sigma_8,NULL,1e-4,1e6,NULL,2048);
	  res = 4.5/pi_sqr*integral;   //see PD97, eq. 29	//4.5/pi_sqr=3*3(in filter)*4pi/(2pi)^3

	  N_SPEC = cosmopar->n;
	  GAMMA = cosmopar->Gamma;
	}
	assert(res>0.0);
	return res;
}


double P_L(double a, double k)	//k in Hubble length-1	//P_L: dimension 1, in unit of Hubble length
{
	double d, t, s;
	const double coverh=ckms/100.;

	if (k==0.0) return(0.0); // avoid nans

	s = sigma_8_sqr();
	d = growfac(a,1.,1.)/growfac(1.,1.,1.);  //numerically solves DGL    
	if((cosmopar->transfer_EH == 1) && !cosmopar->transfer_tabulated) {
	  t = Tsqr_EH(k);
	}
	else if((cosmopar->transfer_EH == 2) && !cosmopar->transfer_tabulated) {
	  t = Tsqr_EH_wiggle(k);
	}
	else if((cosmopar->transfer_EH == 1) && cosmopar->transfer_tabulated) {
	  t = Tsqr_tabulated(k);
	}
	else if (cosmopar->transfer_EH == 0) {   //BBKS
	  t = Tsqr(k);	//Tsqr: (1e-4,1e3)*3000
	}
	else if (cosmopar->transfer_EH == 3) {                                //EBW
	  double knorm = k/coverh;
	  return(Delta_L_smith(knorm)*(d*d)/(k*k*k)*(2.*PI*PI));//Delta_L_smith: (1e-4, 1e3)
	}
	else {
	  printf("Error in routine P_L: no valid transfer function -> %i.\n",cosmopar->transfer_EH);
	  exit(-1);
	}
	return cosmopar->sigma8*cosmopar->sigma8*d*d*t/s;	//initial ps in Tsqr
}


/* PD 22 */
double n_L(double a, double k)
{
	double diff, hh, err, n;

	hh    = k/20.;
	diff = ps_dfridr(P_L, 0.5*k, hh, &err, a); //derivative, defined in stuff.c
	n = .5*k/P_L(a,.5*k)*diff;
	return n;
}


/* CPT 29, PD 15+16 */
double g_PD(double a)
{
  double omega_m, omega_v, fac, a3;

  a3 = a*a*a;
  fac = a + cosmopar->omm*(1.-a) + cosmopar->omv*(a3-a);
  omega_m = cosmopar->omm/fac;
  omega_v = cosmopar->omv*a3/fac;
  return(2.5*omega_m/(pow(omega_m,4./7.) - omega_v + (1. + 0.5*omega_m)*(1. + omega_v/70.)));
}


/* PD 21, 23-27 */
double f_NL(double x, double a, double k)
{
	double A, B, alpha, beta, V;
	double c, gg, f0, f1, f2, f3, f4;

	c = 1.+n_L(a,k)/3.;
	A = .482/pow(c,.947);
	B = .226/pow(c,1.778);
	alpha = 3.31/pow(c,.244);
	beta = .862/pow(c,.287);
	V = 11.55/pow(c,.423);
	gg = g_PD(a);

	f0 = pow(A*x, alpha);
	f1 = 1. + B*beta*x + pow(f0, beta);
	f2 = f0*gg*gg*gg/(V*sqrt(x));
	f3 = 1. + pow(f2, beta);
	f4 = x*pow(f1/f3, 1./beta);
	return f4;
}
 

double P_NL(double a, double k_NL)	// k_NL=(1e-4,1e3)*3000, in Hubble length-1	//P_NL: dimension 1, in unit of Hubble length
{
  const double coverh=ckms/100.;
  if (cosmopar->psnltype!=3) { // use internal power spectrum calculation
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double N_SPEC  = -42.;
	static double GAMMA   = -42.;
	static double SIGMA_8 = -42.;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table_k[N_k], table_P[N_k], table_klog[N_k];
	static double **table_P_NL = 0;
	static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
	double Delta_NL, Delta_L, k_L, lnk_NL=0.0;

	double	omm=0.0,omv=0.0,amp=0.0,ampsqr;
	double  logr1,logr2,diff,rmid,sig,d1,d2,step=5,logr1start,logr2start;
	double	rknl=0.0,rneff=0.0,rncur=0.0;

	double 	sig1,d11,d21;

	double aa, klog, val,logrmidtmp;
	int i,j,niter,itermax=50,wintflag=0,count=0;
	double DIFFACC=0.00001;  //originally 0.001

	if ((OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || SIGMA_8 != cosmopar->sigma8 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH))
	{
	  if(N_a<=1||N_k<=1) printf("error:wrong input for N_a or N_k\n");
	  if (!table_P_NL) {
	    table_P_NL=malloc(N_a*sizeof(double *));
	    for (i=0;i<N_a;i++) {
	      table_P_NL[i]=calloc(N_k,sizeof(double));
	    }
	  }
	  da = (1. - a_min)/(N_a-1.);
	  aa = a_min;
	  logkmin = log(k_min);		//(k_min, k_max)=(1e-3, 1e8)	
	  logkmax = log(k_max);
	  dk =  (logkmax - logkmin)/(N_k-1.); 

	  //FILE *tf = fopen("smitest","w");
		
	  for (i=0; i<N_a; i++, aa +=da) {
	    klog = logkmin;
	    if(cosmopar->psnltype) {
	      omm=om_m(aa);
	      omv=om_v(aa);
	      amp=growfac(aa,omm,omv)/growfac(1.0,cosmopar->omm,OMV_MACRO(1.));	
	      ampsqr=amp*amp;

	      if(aa>0.1) {
		logr1=-2.0;
		logr2=3.5;
		//printf("Start boundaries: %g ... %g\n",(logr1), (logr2));
		count=0;
	      iterstart:

		logr1start=logr1;
		logr2start=logr2;
		      
		niter=0;

		do {
		  rmid=0.5*(logr2+logr1);
		  rmid=pow(10.,rmid);
		  if(wintflag==0) wint1_romb1(rmid,&sig,ampsqr);
		  else wint(rmid,&sig1,&d11,&d21,ampsqr);
		  niter++;
		  diff=sig-1.0;
		  
		  if(diff>DIFFACC) logr1=dlog(rmid);
					
		  if(diff<-DIFFACC) logr2=dlog(rmid);
		} while(fabs(diff)>=DIFFACC && niter<itermax);
			
		if(fabs(diff)>=DIFFACC && niter>=itermax) {
		  count++;
		  if (count>100) DIFFACC*=10.; 
		  if (DIFFACC>0.001) {
		    printf("Error in routine P_NL: Cannot reach accuracy goal in non-linear k loop!\n");
		    exit(-1);
		  } // added by Benjamin to avoid infinite loops

		  logrmidtmp=(logr2start+logr1start)*0.5;
		  if(diff<-DIFFACC) {
		    logr1=logr1start-step;
		    logr2=logr1start;
		  }
		  else if(diff>DIFFACC) {
		    logr1=logr2start;
		    logr2=logr2start+step;
		  }

		  if(logr1<-70 || logr2>70) {
		    wintflag++;
		    if(wintflag>=2) {
		      printf("Couldn't find cosmopar->nonlinear scale int  halofit!\n"); 
		      exit(-99);
		    }
		    printf("Too many bisections - switching to original wint...\n");
		  }
					
		  goto iterstart;
		}

		//...then eff. spectral qties
		if(wintflag==0) wint_romb(rmid,&sig,&d1,&d2,ampsqr);
		else wint(rmid,&sig1,&d11,&d21,ampsqr);
		rknl=1./rmid;	
		rneff=-3.-d1;
		rncur=-d2;
	      }				
	    }

	    //printf("%g  %g  %g  %g\n",1./aa-1.,rneff,rncur,rknl);
			
	    for (j=0; j<N_k; j++, klog += dk) {
	      k_L      = exp(klog);	

	      if(cosmopar->psnltype==0)	{
		Delta_L  = P_L(aa,k_L)*k_L*k_L*k_L/(2.*pi_sqr);	//k_L (1e-3,1e8)
		Delta_NL = f_NL(Delta_L, aa, k_L);
		/* table_k = ln k_{NL}, k_{NL} according to PD 5 */
		lnk_NL = klog + 1./3.*log(1 + Delta_NL);  //PD ~ stable clustering
	      }
	      else {
		Delta_L=amp*amp*Delta_L_smith(k_L/coverh);
		if(aa<=0.1) Delta_NL=Delta_L;
		if(aa>=0.1) {
		  halofit(k_L/coverh,rneff,rncur,rknl,Delta_L,omm,omv,&Delta_NL,aa);	//halofit: (1e-3,1e8)/3000
		  Delta_NL=Delta_NL;	//?!
		}
		lnk_NL=klog;  //SMITH no scaling wrt k //lnk_NL: (1e-3,1e8)
	      }
	      table_k[j] = lnk_NL;
	      if((j>0)&&(table_k[j]<=table_k[j-1])) {
		printf("Error: Interpolation in P_NL does not have monotonically increasing x values for index %i: %10.7g - %10.7g\n",j,table_k[j-1],table_k[j]);
		exit(-1);
	      }
	      /* table_P = ln P_NL, according to PD 3 */
	      table_P[j] = log(2*pi_sqr*Delta_NL) - 3.*lnk_NL;
	      table_klog[j]=klog;  //new array
	    }
	    interpol_spline(table_k,table_P,N_k,table_klog,table_P_NL[i],N_k);

	    /*
	    ps_spline(table_k-1, table_P-1, N_k, cosmopar->n, -2.5, y2-1);

	    klog = logkmin;
	    for (j=0; j<N_k; j++, klog += dk) {
	      ps_splint(table_k-1, table_P-1, y2-1, N_k, klog, &val);
	      table_P_NL[i][j] = val;
	    }
	    */
	    //if(aa>0.1) fprintf(tf,"%g %g %g %g\n",aa,rneff,rncur,rknl);
	  }
	  //fclose(tf);

	  OMEGA_M = cosmopar->omm;
	  OMEGA_V = cosmopar->omv;
	  W0      = cosmopar->w0;
	  WA      = cosmopar->wa;
	  N_SPEC  = cosmopar->n;
	  GAMMA   = cosmopar->Gamma;
	  SIGMA_8 = cosmopar->sigma8;
	  PSNLTYPE	=	cosmopar->psnltype;
	  EH	=	cosmopar->transfer_EH;
	}
	klog = log(k_NL);

	val = interpol2d(table_P_NL, N_a, a_min, 1., da, a, N_k, logkmin, logkmax, dk, klog,cosmopar->n, -2.5); 
	return exp(val);
  }
  else {   // read external power spectrum from file and interpolate
    double logps;
    static int nz_psfile,nk_psfile;
    static double *z_psfile,*k_psfile,**ps_psfile;

    // read file [1st line: # z1 z2 ...; columns: k PS(z1) PS(z2) ...]
    static int flag=0;
    if (!flag) {
      int i,j;
      char dummy[3];
      FILE *dat;
      nz_psfile=get_file_columns(cosmopar->extpsfile)-1;
      nk_psfile=get_file_rows(cosmopar->extpsfile)-1;

      z_psfile=calloc(nz_psfile,sizeof(double));
      k_psfile=calloc(nk_psfile,sizeof(double));
      ps_psfile=calloc(nz_psfile,sizeof(double *));
      for(j=0;j<nz_psfile;j++) {
	ps_psfile[j]=calloc(nk_psfile,sizeof(double));
      }

      if ((dat=fopen(cosmopar->extpsfile,"r"))==NULL) {
	printf("In routine P_NL: couldn't open file %s\n",cosmopar->extpsfile);
	exit(-1);
      }
      fscanf(dat,"%s",dummy);   // initial '#'
      for(j=0;j<nz_psfile;j++) {
	fscanf(dat,"%lf",&z_psfile[j]);
      }
      for(i=0;i<nk_psfile;i++) {
	fscanf(dat,"%lf",&k_psfile[i]);   
	k_psfile[i]=log(k_psfile[i]*coverh);  // transform to k in units of Hubble length
	for(j=0;j<nz_psfile;j++) {
	  fscanf(dat,"%lf",&ps_psfile[j][i]);  
	  ps_psfile[j][i]=log(ps_psfile[j][i]);
	}
      }
      fclose(dat);
      flag=1;
    }

    // interpolate
    logps=interpol_linear_extra_2D(z_psfile,nz_psfile,k_psfile,nk_psfile,ps_psfile,1./a-1.,log(k_NL));
    return(exp(logps)/ps_ch_rescale);
  }	
}


// hubble parameter
double hubble(double a)  // returns H/c in h/Mpc
{
  double matter=cosmopar->omm/(a*a*a);
  double curvature=(1.-cosmopar->omm-cosmopar->omv)/(a*a);
  double vacuum=OMV_MACRO(a);
  const double coverh=ckms/100.;
  return(1./coverh*sqrt(matter+curvature+vacuum));
}


double int_for_w(double a, void *args)
{
  double asqr;
  asqr = a*a;
  return 1./sqrt(a*cosmopar->omm + asqr*(1.-cosmopar->omm-cosmopar->omv) + asqr*asqr*OMV_MACRO(a));
}


/* BS01 2.41, with a(z_2) = 0, a(z_1) = z, in units of Hubble radius c/H_0 = 3000 / h * Mpc */
double w(double a)	//a -> kai
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0      = -42.;
  static double WA      = -42.;
  static double table[N_a];
  static double da = 0.0;
  double aa;
  int i;

  if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa) {
    da = (1.-a_min)/(N_a-1.);
    aa = a_min;
    for (i=0; i<N_a-1; i++, aa+=da) {
      table[i]=int_GSL_integrate_qag(int_for_w,NULL,aa,1.,NULL,2048);
    }
    table[N_a-1] = .0;
    OMEGA_M = cosmopar->omm;
    OMEGA_V = cosmopar->omv;
    W0      = cosmopar->w0;
    WA      = cosmopar->wa;
  }
  return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}


/* BS01 2.4, 2.30 */
double f_K(double w)
{
  double K, K_h, f;

  K = cosmopar->omm + cosmopar->omv - 1;
  if (K > epsilon3) {           /* open */
    K_h = sqrt(K);
    f = 1./K_h*sin(K_h*w);
  } 
  else if (K < -epsilon3) {     /* closed */
    K_h = sqrt(-K);
    f = 1./K_h*sinh(K_h*w);
  } 
  else {                        /* flat */
    f = w;
  }
  return f;
}



/* S98 2.11 */ //beta_p=0: p(z)=delta(z-bg_zdistrpar->z0)
double prob(double z)	//source gal distrib.	//using distrib. parameters in bg_zdistrpar->xx
//normalized to 1??
{
	static double norm = 0.;
	static double BETA_P = -42.;
	static double ALPHA = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	double x, f;

	//First, compute the normalization
	if (ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax)
	{
		if(bg_zdistrpar->beta>0)
		{
			norm	=	int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,1024);
// 			norm = 1/(ps_qromb1(int_for_prob,bg_zdistrpar->zmin,bg_zdistrpar->zmax));
		}
		else
		{
			printf("  Error: sth went wrong during zdistr normalization\n");
			exit(1);
		}
		BETA_P = bg_zdistrpar->beta;
		ALPHA  = bg_zdistrpar->alpha;
		Z0     = bg_zdistrpar->z0;
		ZMIN   = bg_zdistrpar->zmin;
		ZMAX   = bg_zdistrpar->zmax;
	}
	
	if(bg_zdistrpar->beta==0)
	{
		if(fabs(z-bg_zdistrpar->z0)<1e-2)
			return 1.0;
		else
			return 0.0;
	}
	else
	{
		if((bg_zdistrpar->zmin || bg_zdistrpar->zmax) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin))
			return 0.0;
	
		if (z<-epsilon3)
			error("negative z in prob");
		
		x = z/bg_zdistrpar->z0;
		
		return pow(x,bg_zdistrpar->alpha)*exp(-pow(x,bg_zdistrpar->beta))/norm;
	}

	printf("  Error in prob()\n");
	return 0;
}


double int_for_g(double aprime, void *args)
{
	double ww, wprime;
	ww = w(aglob);	//aglob is to be integrated over
	wprime = w(aprime);
	return prob(1./aprime-1.)*f_K(wprime-ww)/f_K(wprime)/(aprime*aprime);
}


/* S98 2.9 */
double g_source(double a)	//the lensing efficiency of lens at a to the gal sample p(z)~prob(1./aprime-1.)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double BETA_P  = -42.;
	static double ALPHA  = -42.;
	static double Z0      = -42.;
	static double table[N_a];
	static double da = 0.0;
	double   aa, wdelta;
	int    i;
 
	/* single source redshift at bg_zdistrpar->z0, p(w) = delta(w-w0) */
	if (bg_zdistrpar->beta == 0.0) {
		aa = 1./(1.+bg_zdistrpar->z0);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}
	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta|| Z0 != bg_zdistrpar->z0) {
		da = (1.-a_min)/(N_a-1.);
		table[0] = 0.0;
		aa = a_min+da;
		for (i=1;i<N_a-1;i++,aa+=da) {
			aglob = aa;
			table[i] = int_GSL_integrate_qag(int_for_g,NULL,a_min,aa,NULL,2048);//ps_qromb(int_for_g, a_min, aa);
		}
		table[N_a-1] = 1.;
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		BETA_P  = bg_zdistrpar->beta;
		ALPHA   = bg_zdistrpar->alpha;
		Z0      = bg_zdistrpar->z0;
	}
	return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}


/* ! global variable sglob is needed here ! */
double int_for_p_2(double a, void *args)
{
  double hoverh0, asqr, g, s, fKw, f, res;

  if (a >= 1.0) error("a>=1 in int_for_p_2");
  s       = sglob;
  asqr    = a*a;
  fKw     = f_K(w(a));
  f       = s/fKw;

  hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));
  g = g_source(a);
  res = g*g/(asqr*asqr)/hoverh0;
  if (cosmopar->nonlinear) res *= P_NL(a, f);
  else res *= P_L(a, f);
  return res;
}


/* S98  3.4 */
double P_2(double s)	//lensing powerspec. 2d without binning
{//time is wasted: calculating for all s then interpolate 
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double N_SPEC  = -42.;
	static double GAMMA   = -42.;
	static double ALPHA  = -42.;
	static double BETA_P  = -42.;
	static double ZMAX    = -42.;
	static double ZMIN    = -42.;
	static double Z0      = -42.;
	static double NONLINEAR = -42;
	static double SIGMA_8 = -42.;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table[N_s];
	static double ds = .0, logsmin = .0, logsmax = .0;
	double   ss, slog, f1, f2,integral;
	int    i;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta ||  ALPHA != bg_zdistrpar->alpha || Z0 != bg_zdistrpar->z0 ||		NONLINEAR != cosmopar->nonlinear || SIGMA_8 != cosmopar->sigma8 || ZMAX != bg_zdistrpar->zmax || ZMIN != bg_zdistrpar->zmin || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		logsmin = log(s_min);
		logsmax = log(s_max);
		ds = (logsmax - logsmin)/(N_s - 1.);
		slog = logsmin;
		for (i=0; i<N_s; i++, slog+=ds)
		{
			ss = exp(slog);
			sglob = ss;

			integral=int_GSL_integrate_qag(int_for_p_2,NULL,a_min,1.0,NULL,2048);
// 			f1 = ps_qromb(int_for_p_2, a_min, 0.7);
// 			f2 = ps_qromo(int_for_p_2, .7, 1.0, ps_midpnt);
			table[i] = log(9./4.*cosmopar->omm*cosmopar->omm*integral);
		}
		
		OMEGA_M =   cosmopar->omm;
		OMEGA_V =   cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		N_SPEC  =   cosmopar->n;
		GAMMA   =   cosmopar->Gamma;
		BETA_P  =   bg_zdistrpar->beta;
		ALPHA  =   bg_zdistrpar->alpha;
		Z0      =   bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		SIGMA_8 =  cosmopar->sigma8;
		ZMAX    =  bg_zdistrpar->zmax;
		ZMIN    =  bg_zdistrpar->zmin;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
	}
	slog = log(s);
	f1 = interpol(table, N_s, logsmin, logsmax, ds, slog, cosmopar->n, cosmopar->n-4.0);

	return exp(f1);
}


double int_for_prob(double z, void *args)
{
	double zz=z/bg_zdistrpar->z0;
	return pow(zz,bg_zdistrpar->alpha)*exp(-pow(zz,bg_zdistrpar->beta));
}


double gauss(double z,double zc,double sigma)
{
  return(exp(-(z-zc)*(z-zc)/(2.*sigma*sigma))/(sigma*sqrt(2.*PI)));
}


double int_for_prob_photoz_inner(double z, void *args)
{
        double *p=(double *)args;
	double zpl=p[0]+bg_zdistrpar->photo_deltaz;
	double zmi=p[0]-bg_zdistrpar->photo_deltaz; 
	return((1.-bg_zdistrpar->photo_fcat)*gauss(z,p[0],bg_zdistrpar->photo_sig*(1.+p[0]))+bg_zdistrpar->photo_fcat*(gauss(z,zpl,bg_zdistrpar->photo_sig*(1.+zpl))+gauss(z,zmi,bg_zdistrpar->photo_sig*(1.+zmi))));   //+zmi !
}

double int_for_prob_photoz1(double z, void *args)
{
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static double ALPHA  = -42.;
	static double BETA_P = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	static int BIN       = -42.;
	static double table[N_z];
	static double dz=0.0;
	static int FLAG = 0;
	static double norm[N_z];

	int i;
	double zz;
        double *p=(double *)args;
	double par[1];

	//gsl_set_error_handler_off();  

	if (FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax  || BIN != (int)p[0])
	{
	  dz = (bg_zdistrpar->zmax-bg_zdistrpar->zmin)/(N_z-1.);
	  zz = bg_zdistrpar->zmin;

  	  for (i=0;i<N_z;i++) {
	    if (zz<1.e-4) table[i]=0.0;
	    else {
	      par[0]=zz;
	      table[i]=int_GSL_integrate_qag(int_for_prob_photoz_inner,par,bg_zdistrpar->tomobin[(int)p[0]],bg_zdistrpar->tomobin[(int)p[0]+1],NULL,2048)*int_for_prob(zz,NULL);		
          if (!FLAG) norm[i]=1./int_GSL_integrate_qag(int_for_prob_photoz_inner,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);  //normalisation for every zz
		  table[i]*=norm[i];
		}
	    zz+=dz;
	  }
	  FCAT   = bg_zdistrpar->photo_fcat;
	  SIG    = bg_zdistrpar->photo_sig;
	  DELTAZ = bg_zdistrpar->photo_deltaz;
  	  ALPHA  = bg_zdistrpar->alpha;
	  BETA_P = bg_zdistrpar->beta;
	  Z0     = bg_zdistrpar->z0;
	  ZMIN   = bg_zdistrpar->zmin;
	  ZMAX   = bg_zdistrpar->zmax;
	  BIN    = (int)p[0];
	  FLAG	 = 1;
	}
	return(interpol(table,N_z,bg_zdistrpar->zmin,bg_zdistrpar->zmax,dz,z,0.,0.));
}

double int_for_prob_photoz2(double z, void *args)
{
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static double ALPHA  = -42.;
	static double BETA_P = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	static int BIN       = -42.;
	static double table[N_z];
	static double dz=0.0;
	static double FLAG = 0;
	static double norm[N_z];

	int i;
	double zz;
        double *p=(double *)args;
	double par[1];

	//gsl_set_error_handler_off();  

	if (FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax  || BIN != (int)p[0])
	{
	  dz = (bg_zdistrpar->zmax-bg_zdistrpar->zmin)/(N_z-1.);
	  zz = bg_zdistrpar->zmin;

  	  for (i=0;i<N_z;i++) {
	    if (zz<1.e-4) table[i]=0.0;
	    else {
	      par[0]=zz;
	      table[i]=int_GSL_integrate_qag(int_for_prob_photoz_inner,par,bg_zdistrpar->tomobin[(int)p[0]],bg_zdistrpar->tomobin[(int)p[0]+1],NULL,2048)*int_for_prob(zz,NULL);		
          if (!FLAG) norm[i]=1./int_GSL_integrate_qag(int_for_prob_photoz_inner,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);  //normalisation for every zz
		  table[i]*=norm[i];
		}
	    zz+=dz;
	  }
	  FCAT   = bg_zdistrpar->photo_fcat;
	  SIG    = bg_zdistrpar->photo_sig;
	  DELTAZ = bg_zdistrpar->photo_deltaz;
  	  ALPHA  = bg_zdistrpar->alpha;
	  BETA_P = bg_zdistrpar->beta;
	  Z0     = bg_zdistrpar->z0;
	  ZMIN   = bg_zdistrpar->zmin;
	  ZMAX   = bg_zdistrpar->zmax;
	  BIN    = (int)p[0];
	  FLAG	 = 1;
	}
	return(interpol(table,N_z,bg_zdistrpar->zmin,bg_zdistrpar->zmax,dz,z,0.,0.));
}

double int_for_prob_photoz3(double z, void *args)
{
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static double ALPHA  = -42.;
	static double BETA_P = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	static int BIN       = -42.;
	static double table[N_z];
	static double dz=0.0;
	static int FLAG = 0;
	static double norm[N_z];

	int i;
	double zz;
        double *p=(double *)args;
	double par[1];

	if (FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax  || BIN != (int)p[0])
	{
	  dz = (bg_zdistrpar->zmax-bg_zdistrpar->zmin)/(N_z-1.);
	  zz = bg_zdistrpar->zmin;

  	  for (i=0;i<N_z;i++) {
	    if (zz<1.e-4) table[i]=0.0;
	    else {
	      par[0]=zz;
	      table[i]=int_GSL_integrate_qag(int_for_prob_photoz_inner,par,bg_zdistrpar->tomobin[(int)p[0]],bg_zdistrpar->tomobin[(int)p[0]+1],NULL,2048)*int_for_prob(zz,NULL);		
          if (!FLAG) norm[i]=1./int_GSL_integrate_qag(int_for_prob_photoz_inner,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);  //normalisation for every zz
		  table[i]*=norm[i];
		}
	    zz+=dz;
	  }
	  FCAT   = bg_zdistrpar->photo_fcat;
	  SIG    = bg_zdistrpar->photo_sig;
	  DELTAZ = bg_zdistrpar->photo_deltaz;
  	  ALPHA  = bg_zdistrpar->alpha;
	  BETA_P = bg_zdistrpar->beta;
	  Z0     = bg_zdistrpar->z0;
	  ZMIN   = bg_zdistrpar->zmin;
	  ZMAX   = bg_zdistrpar->zmax;
	  BIN    = (int)p[0];
	  FLAG	 = 1;
	}
	return(interpol(table,N_z,bg_zdistrpar->zmin,bg_zdistrpar->zmax,dz,z,0.,0.));
}


double probx1(double z,int bin)
{
	static double norm = 0.;
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static double ALPHA  = -42.;
	static double BETA_P = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	static int BIN = -42;
	double x, f;
	double par[1]={bin+0.1};  

	if (bg_zdistrpar->file) {
	      if ((z>bg_zdistrpar->zmin)&&(z<bg_zdistrpar->zmax)) return(probzarb(z,bin));
	      else return(0.0);
	 }

	if (FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax || BIN !=bin)
	{
		if(bg_zdistrpar->beta)
		{
			if(bg_zdistrpar->zmin || bg_zdistrpar->zmax) //Normierung mit Cutoff
			{
			  if(!bg_zdistrpar->photo_sig) //ohne photo-z Fehler
			  {
				norm = 1./int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->tomobin[bin],bg_zdistrpar->tomobin[bin+1],NULL,1024);
			  }
			  else //mit photo-z Fehler
			  {
			        norm = 1./int_GSL_integrate_qag(int_for_prob_photoz1,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);
			  }
			}
			else //Normierung ohne Cutoff
			{
				norm = bg_zdistrpar->beta/(bg_zdistrpar->z0*exp(ps_gammln(3./bg_zdistrpar->beta)));
			}
		}
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		ALPHA  = bg_zdistrpar->alpha;
		BETA_P = bg_zdistrpar->beta;
		Z0     = bg_zdistrpar->z0;
		ZMIN   = bg_zdistrpar->zmin;
		ZMAX   = bg_zdistrpar->zmax;
		BIN = bin;
	}

	if(!bg_zdistrpar->beta)
	{
	        //printf("beta_p=0\n");
		if(fabs(z-bg_zdistrpar->z0)<1e-2) return 1.0;
		else return 0.0;
	}
	if((bg_zdistrpar->photo_sig) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin)) return 0.0;
	if (z<0 && z<-epsilon3) error("negative z in prob");
	if (!bg_zdistrpar->photo_sig)
	{
	  if((bg_zdistrpar->zmin || bg_zdistrpar->zmax) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin || z<bg_zdistrpar->tomobin[bin] || z>bg_zdistrpar->tomobin[bin+1])) return 0.0;
	  if (z<0) x = 0;
	  else x = z/bg_zdistrpar->z0;
	  if (fabs(bg_zdistrpar->beta - 1.) < epsilon3) {
		f = x;
	  } else {
		f = pow(x,bg_zdistrpar->beta);
	  }
	  f=exp(-f);
	  // return norm*DSQR(x)*f;
	  return norm*pow(x,bg_zdistrpar->alpha)*f;
	}
	else
        {
	  return(int_for_prob_photoz1(z,par)*norm);
	}
}


double probx2(double z,int bin)
{
	static double norm = 0.;
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static double ALPHA  = -42.;
	static double BETA_P = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	static int BIN = -42;
	double x, f;
	double par[1]={bin+0.1};  

    if (bg_zdistrpar->file) {
      if ((z>bg_zdistrpar->zmin)&&(z<bg_zdistrpar->zmax)) return(probzarb(z,bin));
      else return(0.0);
    }

	if (FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax || BIN !=bin)
	{
		if(bg_zdistrpar->beta)
		{
			if(bg_zdistrpar->zmin || bg_zdistrpar->zmax) //Normierung mit Cutoff
			{
			  if(!bg_zdistrpar->photo_sig) //ohne photo-z Fehler
			  {
				norm = 1./int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->tomobin[bin],bg_zdistrpar->tomobin[bin+1],NULL,1024);
			  }
			  else //mit photo-z Fehler
			  {
			        norm = 1./int_GSL_integrate_qag(int_for_prob_photoz2,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);
			  }
			}
			else //Normierung ohne Cutoff
			{
				norm = bg_zdistrpar->beta/(bg_zdistrpar->z0*exp(ps_gammln(3./bg_zdistrpar->beta)));
			}
		}
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		ALPHA  = bg_zdistrpar->alpha;
		BETA_P = bg_zdistrpar->beta;
		Z0     = bg_zdistrpar->z0;
		ZMIN   = bg_zdistrpar->zmin;
		ZMAX   = bg_zdistrpar->zmax;
		BIN = bin;
	}
	if(!bg_zdistrpar->beta)
	{
	        //printf("beta_p=0\n");
		if(fabs(z-bg_zdistrpar->z0)<1e-2) return 1.0;
		else return 0.0;
	}
	if((bg_zdistrpar->photo_sig) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin)) return 0.0;
	if (z<0 && z<-epsilon3) error("negative z in prob");
	if (!bg_zdistrpar->photo_sig)
	{
	  if((bg_zdistrpar->zmin || bg_zdistrpar->zmax) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin || z<bg_zdistrpar->tomobin[bin] || z>bg_zdistrpar->tomobin[bin+1])) return 0.0;
	  if (z<0) x = 0;
	  else x = z/bg_zdistrpar->z0;
	  if (fabs(bg_zdistrpar->beta - 1.) < epsilon3) {
		f = x;
	  } else {
		f = pow(x,bg_zdistrpar->beta);
	  }
	  f=exp(-f);
	  //return norm*DSQR(x)*f;
	  return norm*pow(x,bg_zdistrpar->alpha)*f;
	}
	else
        {
	  return(int_for_prob_photoz2(z,par)*norm);
	}
}

double probx3(double z,int bin)
{
	static double norm = 0.;
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static double ALPHA  = -42.;
	static double BETA_P = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	static int BIN = -42;
	double x, f;
	double par[1]={bin+0.1};  

    if (bg_zdistrpar->file) {
      if ((z>bg_zdistrpar->zmin)&&(z<bg_zdistrpar->zmax)) return(probzarb(z,bin));
      else return(0.0);
    }

	if (FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax || BIN !=bin)
	{
		if(bg_zdistrpar->beta)
		{
			if(bg_zdistrpar->zmin || bg_zdistrpar->zmax) //Normierung mit Cutoff
			{
			  if(!bg_zdistrpar->photo_sig) //ohne photo-z Fehler
			  {
				norm = 1./int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->tomobin[bin],bg_zdistrpar->tomobin[bin+1],NULL,1024);
			  }
			  else //mit photo-z Fehler
			  {
			        norm = 1./int_GSL_integrate_qag(int_for_prob_photoz3,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);
			  }
			}
			else //Normierung ohne Cutoff
			{
				norm = bg_zdistrpar->beta/(bg_zdistrpar->z0*exp(ps_gammln(3./bg_zdistrpar->beta)));
			}
		}
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		ALPHA  = bg_zdistrpar->alpha;
		BETA_P = bg_zdistrpar->beta;
		Z0     = bg_zdistrpar->z0;
		ZMIN   = bg_zdistrpar->zmin;
		ZMAX   = bg_zdistrpar->zmax;
		BIN = bin;
	}
	if(!bg_zdistrpar->beta)
	{
	        //printf("beta_p=0\n");
		if(fabs(z-bg_zdistrpar->z0)<1e-2) return 1.0;
		else return 0.0;
	}
	if((bg_zdistrpar->photo_sig) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin)) return 0.0;
	if (z<0 && z<-epsilon3) error("negative z in prob");
	if (!bg_zdistrpar->photo_sig)
	{
	  if((bg_zdistrpar->zmin || bg_zdistrpar->zmax) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin || z<bg_zdistrpar->tomobin[bin] || z>bg_zdistrpar->tomobin[bin+1])) return 0.0;
	  if (z<0) x = 0;
	  else x = z/bg_zdistrpar->z0;
	  if (fabs(bg_zdistrpar->beta - 1.) < epsilon3) {
		f = x;
	  } else {
		f = pow(x,bg_zdistrpar->beta);
	  }
	  f=exp(-f);
	  //return norm*DSQR(x)*f;
	  return norm*pow(x,bg_zdistrpar->alpha)*f;
	}
	else
        {
	  return(int_for_prob_photoz3(z,par)*norm);
	}
}

double probx4(double z,int bin)
{
	static double norm = 0.;
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static double ALPHA  = -42.;
	static double BETA_P = -42.;
	static double Z0     = -42.;
	static double ZMIN   = -42.;
	static double ZMAX   = -42.;
	static int BIN = -42;
	double x, f;
	double par[1]={bin+0.1};  

    if (bg_zdistrpar->file) {
      if ((z>bg_zdistrpar->zmin)&&(z<bg_zdistrpar->zmax)) return(probzarb(z,bin));
      else return(0.0);
    }

	if (FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || ALPHA != bg_zdistrpar->alpha || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN !=bg_zdistrpar->zmin || ZMAX !=bg_zdistrpar->zmax || BIN !=bin)
	{
		if(bg_zdistrpar->beta)
		{
			if(bg_zdistrpar->zmin || bg_zdistrpar->zmax) //Normierung mit Cutoff
			{
			  if(!bg_zdistrpar->photo_sig) //ohne photo-z Fehler
			  {
				norm = 1./int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->tomobin[bin],bg_zdistrpar->tomobin[bin+1],NULL,1024);
			  }
			  else //mit photo-z Fehler
			  {
			        norm = 1./int_GSL_integrate_qag(int_for_prob_photoz3,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);
			  }
			}
			else //Normierung ohne Cutoff
			{
				norm = bg_zdistrpar->beta/(bg_zdistrpar->z0*exp(ps_gammln(3./bg_zdistrpar->beta)));
			}
		}
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		ALPHA  = bg_zdistrpar->alpha;
		BETA_P = bg_zdistrpar->beta;
		Z0     = bg_zdistrpar->z0;
		ZMIN   = bg_zdistrpar->zmin;
		ZMAX   = bg_zdistrpar->zmax;
		BIN = bin;
	}
	if(!bg_zdistrpar->beta)
	{
	        //printf("beta_p=0\n");
		if(fabs(z-bg_zdistrpar->z0)<1e-2) return 1.0;
		else return 0.0;
	}
	if((bg_zdistrpar->photo_sig) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin)) return 0.0;
	if (z<0 && z<-epsilon3) error("negative z in prob");
	if (!bg_zdistrpar->photo_sig)
	{
	  if((bg_zdistrpar->zmin || bg_zdistrpar->zmax) && (z>bg_zdistrpar->zmax || z<bg_zdistrpar->zmin || z<bg_zdistrpar->tomobin[bin] || z>bg_zdistrpar->tomobin[bin+1])) return 0.0;
	  if (z<0) x = 0;
	  else x = z/bg_zdistrpar->z0;
	  if (fabs(bg_zdistrpar->beta - 1.) < epsilon3) {
		f = x;
	  } else {
		f = pow(x,bg_zdistrpar->beta);
	  }
	  f=exp(-f);
	  //return norm*DSQR(x)*f;
	  return norm*pow(x,bg_zdistrpar->alpha)*f;
	}
	else
        {
	  return(int_for_prob_photoz3(z,par)*norm);
	}
}



double int_for_gx1(double aprime, void *args)
{
	double ww, wprime;
	int bin;
	double	*arg	=	(double*)args;

	bin=(int) arg[0];

	ww = w(aglob);
	wprime = w(aprime);
	return probx1(1./aprime-1.,bin)*f_K(wprime-ww)/f_K(wprime)/(aprime*aprime);
}

double int_for_gx2(double aprime, void *args)
{
	double ww, wprime;
	int bin;
	double	*arg	=	(double*)args;

	bin=(int) arg[0];

	ww = w(aglob);
	wprime = w(aprime);
	return probx2(1./aprime-1.,bin)*f_K(wprime-ww)/f_K(wprime)/(aprime*aprime);
}

double int_for_gx3(double aprime, void *args)
{
	double ww, wprime;
	int bin;
	double	*arg	=	(double*)args;

	bin=(int) arg[0];

	ww = w(aglob);
	wprime = w(aprime);
	return probx3(1./aprime-1.,bin)*f_K(wprime-ww)/f_K(wprime)/(aprime*aprime);
}

double int_for_gx4(double aprime, void *args)
{
	double ww, wprime;
	int bin;
	double	*arg	=	(double*)args;

	bin=(int) arg[0];

	ww = w(aglob);
	wprime = w(aprime);
	return probx4(1./aprime-1.,bin)*f_K(wprime-ww)/f_K(wprime)/(aprime*aprime);
}


double g_sourcex1(double a, int bin)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
        static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static int BIN =-42;
	static double table[N_a];
	static double da = 0.0;
	double  arg[2];
	double   aa, wdelta;
	int    i;

	// single source redshift at median redshifts of bins
	if (bg_zdistrpar->sheet) {
		aa = 1./(1.+bg_zdistrpar->zmed[bin]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	// single source redshift at bg_zdistrpar->tomobin[1]
	if (!bg_zdistrpar->beta) {
		aa = 1./(1.+bg_zdistrpar->tomobin[1]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || BIN != bin)
	{
		da = (1.-a_min)/(N_a-1.);
		table[0] = 0.0;
		aa = a_min+da;
		arg[0]=bin+0.1;

		for (i=1;i<N_a-1;i++,aa+=da)
		{
			long double sum=0,step=0.02,limit;

			aglob = aa;
			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)) limit = 1./(bg_zdistrpar->tomobin[bin+1]+1.);
			//else limit = a_min;
			else limit = 1./(1.+bg_zdistrpar->zmax);

			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)&&(aa<1./(bg_zdistrpar->tomobin[bin+1]+1.))) table[i]=0;
			else {
			  while(limit  < aa-step)//!!
			  {
				sum+=int_GSL_integrate_qag(int_for_gx1,arg,limit,limit+step,NULL,1024);
				limit+=step;
			  }
			  sum+=int_GSL_integrate_qag(int_for_gx1,arg,limit,aa,NULL,1024);
			  table[i]=sum;
			}
		}

		table[N_a-1] = 1.;
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		BIN = bin;
	}
	return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}


double g_sourcex2(double a, int bin)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
    static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static int BIN =-42;
	static double table[N_a];
	static double da = 0.0;
	double  arg[2];
	double   aa, wdelta;
	int    i;

	// single source redshift at median redshifts of bins
	if (bg_zdistrpar->sheet) {
		aa = 1./(1.+bg_zdistrpar->zmed[bin]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	// single source redshift at bg_zdistrpar->tomobin[1]
	if (!bg_zdistrpar->beta) {
		aa = 1./(1.+bg_zdistrpar->tomobin[1]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || BIN != bin)
	{
		if(N_a==1) {printf("caution: N_a=1!\n");da = 0.;} else
		da = (1.-a_min)/(N_a-1.);
		table[0] = 0.0;
		aa = a_min+da;
		arg[0]=bin+0.1;

		for (i=1;i<N_a-1;i++,aa+=da)	//for g(a)
		{
			long double sum=0,step=0.02,limit;
			
			aglob = aa;
			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)) limit = 1./(bg_zdistrpar->tomobin[bin+1]+1.);
			//else limit = a_min;
			else limit = 1./(1.+bg_zdistrpar->zmax);

			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)&&(aa<1./(bg_zdistrpar->tomobin[bin+1]+1.))) table[i]=0;
			else {
  			  while(limit  < aa-step)
			  {
				sum+=int_GSL_integrate_qag(int_for_gx2,arg,limit,limit+step,NULL,1024);
				limit+=step;
			  }
			  sum+=int_GSL_integrate_qag(int_for_gx2,arg,limit,aa,NULL,1024);
			  table[i]=sum;
			}
		}
		
		table[N_a-1] = 1.;
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		BIN = bin;
	}
	return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}

double g_sourcex3(double a, int bin)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
    static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static int BIN =-42;
	static double table[N_a];
	static double da = 0.0;
	double  arg[2];
	double   aa, wdelta;
	int    i;

	// single source redshift at median redshifts of bins
	if (bg_zdistrpar->sheet) {
		aa = 1./(1.+bg_zdistrpar->zmed[bin]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	// single source redshift at bg_zdistrpar->tomobin[1]
	if (!bg_zdistrpar->beta) {
		aa = 1./(1.+bg_zdistrpar->tomobin[1]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || BIN != bin)
	{
		da = (1.-a_min)/(N_a-1.);
		table[0] = 0.0;
		aa = a_min+da;
		arg[0]=bin+0.1;

		for (i=1;i<N_a-1;i++,aa+=da)	//for g(a)
		{
			long double sum=0,step=0.02,limit;
			
			aglob = aa;
			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)) limit = 1./(bg_zdistrpar->tomobin[bin+1]+1.);
			//else limit = a_min;
			else limit = 1./(1.+bg_zdistrpar->zmax);
			
			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)&&(aa<1./(bg_zdistrpar->tomobin[bin+1]+1.))) table[i]=0;
			else {
  			  while(limit  < aa-step)
			  {
				sum+=int_GSL_integrate_qag(int_for_gx3,arg,limit,limit+step,NULL,1024);
				limit+=step;
			  }
			  sum+=int_GSL_integrate_qag(int_for_gx3,arg,limit,aa,NULL,1024);
			  table[i]=sum;
			}
		}
		
		table[N_a-1] = 1.;
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		BIN = bin;
	}
	return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}


double g_sourcex4(double a, int bin)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
    static double FCAT   = -42.;
	static double SIG    = -42.;
	static double DELTAZ = -42.;
	static int BIN =-42;
	static double table[N_a];
	static double da = 0.0;
	double  arg[2];
	double   aa, wdelta;
	int    i;

	// single source redshift at median redshifts of bins
	if (bg_zdistrpar->sheet) {
		aa = 1./(1.+bg_zdistrpar->zmed[bin]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	// single source redshift at bg_zdistrpar->tomobin[1]
	if (!bg_zdistrpar->beta) {
		aa = 1./(1.+bg_zdistrpar->tomobin[1]);
		if (a<=aa) return 0.0;
		wdelta = w(aa);
		return f_K(wdelta-w(a))/f_K(wdelta);
	}

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || BIN != bin)
	{
		da = (1.-a_min)/(N_a-1.);
		table[0] = 0.0;
		aa = a_min+da;
		arg[0]=bin+0.1;

		for (i=1;i<N_a-1;i++,aa+=da)	//for g(a)
		{
			long double sum=0,step=0.02,limit;
			
			aglob = aa;
			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)) limit = 1./(bg_zdistrpar->tomobin[bin+1]+1.);
			//else limit = a_min;
			else limit = 1./(1.+bg_zdistrpar->zmax);
			
			if((!bg_zdistrpar->file)&&(!bg_zdistrpar->photo_sig)&&(aa<1./(bg_zdistrpar->tomobin[bin+1]+1.))) table[i]=0;
			else {
  			  while(limit  < aa-step)
			  {
				sum+=int_GSL_integrate_qag(int_for_gx4,arg,limit,limit+step,NULL,1024);
				limit+=step;
			  }
			  sum+=int_GSL_integrate_qag(int_for_gx4,arg,limit,aa,NULL,1024);
			  table[i]=sum;
			}
		}
		
		table[N_a-1] = 1.;
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		FCAT   = bg_zdistrpar->photo_fcat;
		SIG    = bg_zdistrpar->photo_sig;
		DELTAZ = bg_zdistrpar->photo_deltaz;
		BIN = bin;
	}
	return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}

double int_for_p_2x(double a, void *args)
{
	double hoverh0, asqr, s, fKw, f, res;
	int bin1,bin2;

	double *arg	=	(double*)args;
		
	bin1=(int)arg[0];
	bin2=(int)arg[1];

	if (a >= 1.0) error("a>=1 in int_for_p_2x");
	s       = sglob;
	asqr    = a*a;
	fKw     = f_K(w(a));
	f       = s/fKw;

	hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));
	res = g_sourcex1(a,bin1)*g_sourcex2(a,bin2)/(asqr*asqr)/hoverh0;
	if (cosmopar->nonlinear) res *= P_NL(a, f);
	else res *= P_L(a, f);
	return res;
}


double P_2x(double s, int bin1, int bin2)	//lensing tomo P_2
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double N_SPEC  = -42.;
	static double GAMMA   = -42.;
	static double BETA_P  = -42.;
	static double ZMAX    = -42.;
	static double ZMIN    = -42.;
	static double Z0      = -42.;
	static double FCAT    = -42.;
	static double SIG     = -42.;
	static double DELTAZ  = -42.;
	static double NONLINEAR = -42.;
	static double SIGMA_8 = -42.;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table[N_s];
	static double ds = .0, logsmin = .0, logsmax = .0;
	double   ss, slog, f1, f2;
	int    i;
	double arg[2], integral;
   

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || NONLINEAR != cosmopar->nonlinear || SIGMA_8 != cosmopar->sigma8 || ZMAX != bg_zdistrpar->zmax || ZMIN != bg_zdistrpar->zmin || NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)	
	{
		logsmin = log(s_min);
		logsmax = log(s_max);
		if(N_s==1) ds=0.; else
		ds = (logsmax - logsmin)/(1.*N_s - 1.);
		slog = logsmin;
		arg[0]=bin1+0.1;
		arg[1]=bin2+0.1;

		for (i=0; i<N_s; i++, slog+=ds) {
			ss = exp(slog);
			sglob = ss;
  			integral=int_GSL_integrate_qag(int_for_p_2x,arg,a_min,1.0,NULL,2048);
			table[i] = log(9./4.*cosmopar->omm*cosmopar->omm*integral);
		}
		OMEGA_M =   cosmopar->omm;
		OMEGA_V =   cosmopar->omv;
		W0      =   cosmopar->w0;
		WA      =   cosmopar->wa;
		N_SPEC  =   cosmopar->n;
		GAMMA   =   cosmopar->Gamma;
		BETA_P  =   bg_zdistrpar->beta;
		Z0      =   bg_zdistrpar->z0;
		FCAT    =   bg_zdistrpar->photo_fcat;
		SIG     =   bg_zdistrpar->photo_sig;
		DELTAZ  =   bg_zdistrpar->photo_deltaz;
		NONLINEAR = cosmopar->nonlinear;
		SIGMA_8 =   cosmopar->sigma8;
		ZMAX    =   bg_zdistrpar->zmax;
		ZMIN    =   bg_zdistrpar->zmin;
		NBIN1 = bin1;
		NBIN2 = bin2;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
	}
	slog = log(s);
	f1 = interpol(table, N_s, logsmin, logsmax, ds, slog, cosmopar->n, cosmopar->n-4.0);
	return exp(f1);
}


// identical to P_2x, but uses external 3D matter power spectrum
double int_for_p_2xext(double a, void *args)
{
  double hoverh0, asqr, s, fKw, f, res, g1, g2;
  int bin1,bin2,na,nk;
  const double coverh=ckms/100.;

  double *arg=(double*)args;
  bin1=(int)arg[0];
  bin2=(int)arg[1];
  na=(int)arg[2];
  nk=(int)arg[3];

  static double *aa,*kk;
  static double **ps;

  static int flag=0;
  if (!flag) {
    int i,j;
    double dummy;
    FILE *dat;
    aa=calloc(na,sizeof(double));
    kk=calloc(nk,sizeof(double));
    ps=calloc(na,sizeof(double *));
    for(j=0;j<na;j++) {
      ps[j]=calloc(nk,sizeof(double));
    }

    if ((dat=fopen(extpsfile,"r"))==NULL) {
      printf("Couldn't open file %s\n",extpsfile);
      exit(-1);
    }
    fscanf(dat,"%lf",&dummy);
    for(j=0;j<na;j++) {
      fscanf(dat,"%lf",&aa[j]);
    }
    for(i=0;i<nk;i++) {
      fscanf(dat,"%lf",&kk[i]);   // log_10(k)
      kk[i]+=log10(coverh);  // transform to k in units of Hubble length
      for(j=0;j<na;j++) {
	fscanf(dat,"%lf",&ps[j][i]);  //log_10(P)
      }
    }
    fclose(dat);
    flag=1;
  }

  if (a >= 1.0) error("a>=1 in int_for_p_2xext");
  s       = sglob;
  asqr    = a*a;
  fKw     = f_K(w(a));
  f       = s/fKw;   

  if (bg_zdistrpar->sheet) {
    double asource,wdelta;
    asource=1./(1.+bg_zdistrpar->tomobin[bin1+1]);
    if (a<=asource) g1=0.0;
    else {
      wdelta=w(asource);
      g1=f_K(wdelta-w(a))/f_K(wdelta);
    }

    asource=1./(1.+bg_zdistrpar->tomobin[bin2+1]);
    if (a<=asource) g2=0.0;
    else {
      wdelta=w(asource);
      g2=f_K(wdelta-w(a))/f_K(wdelta);
    }
  }
  else {
    g1=g_sourcex1(a,bin1);
    g2=g_sourcex2(a,bin2);
  }

  if ((g1==0.0)||(g2==0.0)) res=0.0;
  else {
    hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));
    res=g1*g2/(asqr*asqr)/hoverh0*pow(10.,interpol_linear_ext(aa,na,kk,nk,ps,a,log10(f)))/ps_ch_rescale;   //convert 3dps to units of Hubble length
  }
  return res;
}


double P_2xext(double s, int bin1, int bin2, int zbins, int kbins)
{
  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0      = -42.;
  static double WA      = -42.;
  static double N_SPEC  = -42.;
  static double GAMMA   = -42.;
  static double BETA_P  = -42.;
  static double ZMAX    = -42.;
  static double ZMIN    = -42.;
  static double Z0      = -42.;
  static double FCAT    = -42.;
  static double SIG     = -42.;
  static double DELTAZ  = -42.;
  static double NONLINEAR = -42.;
  static double SIGMA_8 = -42.;
  static int NBIN1=-42;
  static int NBIN2=-42;
  static int PSNLTYPE=-42;
  static int EH=-42;
  static double table[N_s];
  static double ds = .0, logsmin = .0, logsmax = .0;
  double   ss, slog, f1;
  int    i;
  double arg[4], integral;

  if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || NONLINEAR != cosmopar->nonlinear || SIGMA_8 != cosmopar->sigma8 || ZMAX != bg_zdistrpar->zmax || ZMIN != bg_zdistrpar->zmin || NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH) {
    logsmin = log(s_min);
    logsmax = log(s_max);
    if(N_s==1) ds=0.; 
    else ds = (logsmax - logsmin)/(1.*N_s - 1.);
    slog = logsmin;
    arg[0]=bin1+0.1;
    arg[1]=bin2+0.1;
    arg[2]=zbins+0.1;
    arg[3]=kbins+0.1;

    for (i=0; i<N_s; i++, slog+=ds) {
      ss = exp(slog);
      sglob = ss;
      integral=int_GSL_integrate_qag(int_for_p_2xext,arg,a_min,1.0,NULL,2048);
      table[i] = log(9./4.*cosmopar->omm*cosmopar->omm*integral);
    }
    OMEGA_M =   cosmopar->omm;
    OMEGA_V =   cosmopar->omv;
    W0      =   cosmopar->w0;
    WA      =   cosmopar->wa;
    N_SPEC  =   cosmopar->n;
    GAMMA   =   cosmopar->Gamma;
    BETA_P  =   bg_zdistrpar->beta;
    Z0      =   bg_zdistrpar->z0;
    FCAT    =   bg_zdistrpar->photo_fcat;
    SIG     =   bg_zdistrpar->photo_sig;
    DELTAZ  =   bg_zdistrpar->photo_deltaz;
    NONLINEAR = cosmopar->nonlinear;
    SIGMA_8 =   cosmopar->sigma8;
    ZMAX    =   bg_zdistrpar->zmax;
    ZMIN    =   bg_zdistrpar->zmin;
    NBIN1 = bin1;
    NBIN2 = bin2;
    PSNLTYPE	=	cosmopar->psnltype;
    EH	=	cosmopar->transfer_EH;
  }
  slog = log(s);
  f1 = interpol(table, N_s, logsmin, logsmax, ds, slog, cosmopar->n, cosmopar->n-4.0);
  return exp(f1);
}





/*************************************************************************************/
/*  Functions for AIFA-HALOFIT                                                                                                     */
/*************************************************************************************/

double Delta_L_smith(double k)
{
  const double coverh=ckms/100.;
  double kkorr=k*coverh;
  double k3=kkorr*kkorr*kkorr;
  double fac=cosmopar->sigma8*cosmopar->sigma8*k3/(2.*PI*PI*sigma_8_sqr());
  if (cosmopar->transfer_EH==0) return(fac*Tsqr(kkorr));            //BBKS
  else if (cosmopar->transfer_EH==1) return(fac*Tsqr_EH(kkorr));    //EH
  else if (cosmopar->transfer_EH==2) return(fac*Tsqr_EH_wiggle(kkorr));    //EH w/ wiggles
  else if (cosmopar->transfer_EH==3) {      //original EBW
	double 		keff,q,q8,tk,tk8;
 
	keff=0.172+0.011*log(cosmopar->Gamma/0.36)*log(cosmopar->Gamma/0.36);
	q=1e-20+k/cosmopar->Gamma;
	q8=1e-20+keff/cosmopar->Gamma;
	tk=1/pow(1+pow(6.4*q+pow(3.0*q,1.5)+(1.7*q)*(1.7*q),1.13),(1/1.13));
	tk8=1/pow(1+pow(6.4*q8+pow(3.0*q8,1.5)+(1.7*q8)*(1.7*q8),1.13),(1/1.13));
 
	return cosmopar->sigma8*cosmopar->sigma8*pow(q/q8,3+cosmopar->n)*tk*tk/tk8/tk8;
  }
  else {
    printf("Error in routine Delta_L_smith: no valid transfer function.\n");
    exit(-1);
  }
}


double Delta_L_smith_tab(double k)
{
	double		keff,q,q8,tk,tk8,kmin=1e-6,kmax=1e6,u,v,logkk,logk;
	int			ntab=1000,i;
	static int		flag=0;
	static double	*table,logmin,logmax,dlogk;
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double N_SPEC  = -42.;
	static double GAMMA   = -42.;
	static double SIGMA_8 = -42.;
	static double BETA_P   = -42.;
	static double Z0 = -42.;
	static double NONLINEAR = -42.;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;

		if(!flag)
		{
			table=malloc(ntab*sizeof(double));
			flag=1;
		}

		logmin=log(kmin); logmax=log(kmax);
		dlogk=(logmax-logmin)/(1.0*ntab);
  
		for(i=0;i<ntab;i++)
		{
			logkk=logmin+i*dlogk;
			table[i]=log(Delta_L_smith(exp(logkk)));
		}
	}

	logk=log(k);
	i=(int)((logk-logmin)/dlogk);
	u=(logk-logmin)/dlogk-1.0*i;
	v=1.0-u;

	if(i<ntab-1 && i>=0) return exp(v*table[i] + u*table[i+1]);
	else if(i>=ntab-1 || i<0) return(Delta_L_smith(k));
	return -1;
}


double om_m(double a)
{
  return cosmopar->omm/(cosmopar->omm + a*(a*a*OMV_MACRO(a) + (1-cosmopar->omm-cosmopar->omv)));
}
 
double om_v(double a)
{
  return OMV_MACRO(a)*a*a*a/(cosmopar->omm + a*(a*a*OMV_MACRO(a) + (1-cosmopar->omm-cosmopar->omv)));
}

/*
double growfac(double a, double omm, double omv)
{
        double m=1-6.*cosmopar->w0/5.;
        double alpha=3./(5.-cosmopar->w0/(1-cosmopar->w0))+3./125.*(1-cosmopar->w0)*(1-3.*cosmopar->w0/2.)/(m*m*m)*(1-omm);
	return a*2.5*omm/(pow(omm,alpha)-omv+(1+0.5*omm)*(1+(-0.28/(cosmopar->w0+0.08)-0.3)*omv));   //Percival 2005, only valid for w_a=0 & flat universe!
}
*/



//function for growfac (DGL)
int func_for_growfac(double a,const double y[],double f[],void *params)
{
  double *p=(double *)params;
  if (a == 0) {
    printf("a=0 in function 'func_for_growfac'!\n");
    exit(1);
  }
  double aa=a*a;
  double omegam=p[0]/(aa*a);
  double omegav=p[1]*exp(-3.*((p[2]+p[3]+1)*log(a)+p[3]*(1.-a)));
  double hub=omegam+(1-p[0]-p[1])/(a*a)+omegav;
  f[0]=y[1];
  f[1]=y[0]*3.*p[0]/(2.*hub*aa*aa*a)-y[1]/a*(2.-(omegam+(3.*(p[2]+p[3]*(1.-a))+1)*omegav)/(2.*hub));
  return GSL_SUCCESS;
}

double growfac(double a, double omm, double omv)
{
        //variables omm,omv ignored!
	const double MINA=1.e-8;
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double ai[N_a];
	static double table[N_a];
	double res;

	gsl_interp *intf=gsl_interp_alloc(gsl_interp_linear,N_a);
	gsl_interp_accel *acc=gsl_interp_accel_alloc();

	if (a-1.>0.01) {
	  printf("Error in routine 'growfac': encountered scale factor larger than 1: %g\n",a);
	  exit(-1);
	}

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa)
	{
	  OMEGA_M = cosmopar->omm;
	  OMEGA_V = cosmopar->omv;
	  W0      = cosmopar->w0;
	  WA      = cosmopar->wa;

  	  int i;
	  const gsl_odeiv_step_type *T=gsl_odeiv_step_rkf45;
	  gsl_odeiv_step *s=gsl_odeiv_step_alloc(T,2);
	  gsl_odeiv_control *c=gsl_odeiv_control_y_new(1.e-6,0.0);  //abs./rel. allowed error on y
	  gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(2);

	  double t=MINA;            //start a
	  double t1=1.;                //final a
	  double h=1.e-6;              //initial step size
	  double y[2]={MINA,MINA};   //initial conditions
	  double norm=0.0;

	  double par[4]={OMEGA_M,OMEGA_V,W0,WA};    //Omega_m,Omega_v,w_0,w_a
          gsl_odeiv_system sys={func_for_growfac,NULL,2,&par};

	  for (i=1;i<=N_a;i++) {
	    ai[i-1]=i*t1/(1.*N_a);
	    while(t<ai[i-1]) gsl_odeiv_evolve_apply(e,c,s,&sys,&t,ai[i-1],&h,y);
	    if (i==1) norm=y[0]/ai[i-1];
	    table[i-1]=y[0]/norm;
	  }

	  gsl_odeiv_evolve_free(e);
	  gsl_odeiv_control_free(c);
	  gsl_odeiv_step_free(s);
	}
        gsl_interp_init(intf,ai,table,N_a);
	if (a-1.>0.0) a=1.;  // avoid nans due to round-off errors
	res=gsl_interp_eval(intf,ai,table,a,acc);
	gsl_interp_accel_free(acc);
	gsl_interp_free(intf);
	return(res);
}



#define USE_GAMMA 1       // 2: use Percival approximation; 1: use gamma as parameter; 0: calculate derivative
#define gamma_growth 0.55 // should eventually be implemented as parameter
#define step_a 0.004      // step size for derivative

double growth_rate(double a)
{
  double res=0.0;

  if (USE_GAMMA==2) {
    double omega_m, hub, a3, m, alpha;
    a3 = a*a*a;
    hub=cosmopar->omm/a3+OMV_MACRO(a);
    omega_m = cosmopar->omm/a3/hub;
    m=1-6.*cosmopar->w0/5.;
    alpha=3./(5.-cosmopar->w0/(1-cosmopar->w0))+3./125.*(1-cosmopar->w0)*(1-3.*cosmopar->w0/2.)/(m*m*m)*(1-omega_m);
    res=pow(omega_m,alpha);  //see Percival 2005, only valid for w_a=0 and flat universe!
  }
  else if (USE_GAMMA==1) {
    double omega_m, hub, a3;
    a3 = a*a*a;
    hub=cosmopar->omm/a3+OMV_MACRO(a);
    omega_m = cosmopar->omm/a3/hub;
    res=pow(omega_m,gamma_growth);
  }
  else if (USE_GAMMA==0) {
    if (a+2.*step_a-1.>0.0) {
      printf("Error in routine 'growth_rate': encountered scale factor too close to unity: %g\n",a);
      exit(-1);
    }
    if (a-2.*step_a<0.0) {
      printf("Error in routine 'growth_rate': encountered scale factor too close to 0: %g\n",a);
      exit(-1);
    }

    int i;
    double deriv,D[5];
    for(i=0;i<5;i++) {
      D[i]=growfac(a+(i-2)*step_a,1.,1.);
    }
    deriv=(-D[4]+8.*D[3]-8.*D[1]+D[0])/(12.*step_a);
    res=a/D[2]*deriv;
  }
  else {
    printf ("Error in routine 'growth_rate': incorrect value for USE_GAMMA.\n");
    exit(-1);
  }
  return(res);
}

#undef USE_GAMMA
#undef gamma_growth
#undef step_a


void wint(double r,double *sig,double *d1,double *d2, double ampsqr)
{
	double 		sum1,sum2,sum3,t,y,x,w1,w2,w3,den,xsqr;
	int		i,nint;

	nint=10000;
	sum1=0.0;
	sum2=0.0;
	sum3=0.0;
	for(i=1;i<=nint;i++)
	{
		t=(i-0.5)/(1.0*nint);
		y=-1.0+1.0/t;
		*d2=ampsqr*Delta_L_smith(y);
		x=y*r;
		xsqr=x*x;
		w1=exp(-xsqr);
		w2=2.0*xsqr*w1;
		w3=4.0*xsqr*(1.0-xsqr)*w1;
		den=(*d2)/y/t/t;
		sum1+=w1*den;
		sum2+=w2*den;
		sum3+=w3*den;
	}
	sum1/=(1.0*nint);
	sum2/=(1.0*nint);
	sum3/=(1.0*nint);
	*sig=sqrt(sum1);
	*d1=-sum2/sum1;
	*d2=-sum2*sum2/sum1/sum1 - sum3/sum1;
}

 
 void wint_romb(double r,double *sig,double *d1,double *d2, double ampsqr)
{
	double 		sum1,sum2,sum3,arg[2];
  
	sum1	=	0.0;
	sum2	=	0.0;
	sum3	=	0.0;

	arg[0]=	r;
	arg[1]=	ampsqr;

	sum1	=	int_GSL_integrate_qag(intforwint_knl1,arg,log(1e-4),log(1e5),NULL,1024);
	sum2	=	int_GSL_integrate_qag(intforwint_neff,arg,log(1e-4),log(1e5),NULL,1024);
	sum3	=	int_GSL_integrate_qag(intforwint_cur,arg,log(1e-4),log(1e5),NULL,1024);

	*sig	=	sqrt(sum1);
	*d1	=-	sum2/sum1;
	*d2	=-	sum2*sum2/sum1/sum1 - sum3/sum1;
}


void wint1_romb1(double r, double *sig, double ampsqr)
{
	double		sum,lim,arg[2];

	sum=0;
 
	lim=0.5*(sqrt(5.*log(10.))/r+1.);
	if (lim<1000.) lim=1000.;
	
	arg[0]=r;
	arg[1]=ampsqr;

	sum=int_GSL_integrate_qag(intforwint_knl1,arg,log(1e-4),log(1e5),NULL,1024);

	*sig=sqrt(sum);
}


double intforwint_knl1(double lnk, void *args)
{
 //args: 0: r;; 1: ampsqr
	double k;

	double *arg = (double*) args;
	k=exp(lnk);
	return exp(-k*arg[0]*k*arg[0])*arg[1]*Delta_L_smith_tab(k);
 
}


double intforwint_neff(double lnk, void *args)
{
 //args: 0: r;; 1: ampsqr
	double xsqr,k;

	double *arg = (double*) args;
	k=exp(lnk);
	xsqr=k*arg[0]*k*arg[0];
	return 2.0*xsqr*exp(-k*arg[0]*k*arg[0])*arg[1]*Delta_L_smith_tab(k);
}


double intforwint_cur(double lnk, void *args)
{
 //args: 0: r;; 1: ampsqr
	double xsqr,k;
	double *arg = (double*) args;
	
	k=exp(lnk);
	xsqr=k*arg[0]*k*arg[0];
	return 4.0*xsqr*(1.0-xsqr)*exp(-xsqr)*arg[1]*Delta_L_smith_tab(k);
}


void halofit(double rk,double rn,double rncur,double rknl,double plin,double om_m, double om_v,double *pnl,double aa)   //double aa neu!
{
	double 		gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
	double 		y,w_de;
	double 		f1a,f2a,f3a,f1b,f2b,f3b,frac,pq,ph;


	w_de=cosmopar->w0+cosmopar->wa*(1.-aa);   //w_z(z);

	if (cosmopar->psnltype==1) {  // Smith et al. (2003)
	  a=pow(10,(1.4861+1.8369*rn+1.6762*rn*rn+0.7940*rn*rn*rn+0.1670*rn*rn*rn*rn-0.6206*rncur));
	  b=pow(10,(0.9463+0.9466*rn+0.3084*rn*rn-0.940*rncur));
	  c=pow(10,(-0.2807+0.6669*rn+0.3214*rn*rn-0.0793*rncur));
	  xmu=pow(10,(-3.5442+0.1908*rn));
	  xnu=pow(10,(0.9589+1.2857*rn));
	  alpha=1.38848+0.3700*rn-0.1452*rn*rn;
	  beta=0.8291+0.9854*rn+0.3401*rn*rn;
	  gam=0.86485+0.2989*rn+0.1631*rncur;
	}
	else if (cosmopar->psnltype==2) { // Takahashi et al. (2012)
	  a=pow(10,(1.5222+2.8553*rn+2.3706*rn*rn+0.9903*rn*rn*rn+0.2250*rn*rn*rn*rn-0.6038*rncur+0.1749*om_v*(1.+w_de)));
	  b=pow(10,(-0.5642+0.5864*rn+0.5716*rn*rn-1.5474*rncur+0.2279*om_v*(1.+w_de)));
	  c=pow(10,(0.3698+2.0404*rn+0.8161*rn*rn+0.5869*rncur));
	  xmu=0.0;
	  xnu=pow(10,(5.2105+3.6902*rn));
	  alpha=fabs(6.0835+1.3373*rn-0.1959*rn*rn-5.5274*rncur);
	  beta=2.0379-0.7354*rn+0.3157*rn*rn+1.2490*rn*rn*rn+0.3980*rn*rn*rn*rn-0.1682*rncur;
	  gam=0.1971-0.0843*rn+0.8460*rncur;
	}
	else {
	  printf("Error in routine HALOFIT: Incorrect parameter for non-linear power spectrum correction: %i\n",cosmopar->psnltype);
	  exit(-1);
	}

	if(fabs(1-om_m)>0.01)
	{
		f1a=pow(om_m,(-0.0732));
		f2a=pow(om_m,(-0.1423));
		f3a=pow(om_m,(0.0725));
		f1b=pow(om_m,(-0.0307));
		f2b=pow(om_m,(-0.0585));
		f3b=pow(om_m,(0.0743));
		frac=om_v/(1.-om_m);

		if (cosmopar->psnltype==1) {   //interpolation in eq. of state parameters
		  double we;
		  we=(frac*w_de)+((1.0-frac)*(-1.0/3.0));
		  frac=(-1.0)*((3.0*we)+1.0)/2.0;
		}

		f1=frac*f1b + (1-frac)*f1a;
		f2=frac*f2b + (1-frac)*f2a;
		f3=frac*f3b + (1-frac)*f3a;
	}
	else
	{
		f1=1.0;
		f2=1.;
		f3=1.;
	}
	y=(rk/rknl);
	ph=a*pow(y,f1*3)/(1+b*pow(y,f2)+pow(f3*c*y,3-gam));
	ph=ph/(1+xmu/y+xnu/y/y);
	pq=plin*pow(1+plin,beta)/(1+plin*alpha)*exp(-y/4.0-y*y/8.0);
	*pnl=pq+ph;
}
 
double dlog(double x)
{
	return log(x)/log(10.0);
}




/*************************************************************************************/
/*  Functions for Hankel transform                                                   */
/*************************************************************************************/

//If xi(theta=0) is needed, uncomment line 2275 (massive slowdown)!

double Hankel_Xi(double theta, int pm, char *path)
{
	int              i;

	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double ZMIN 	 	= -42.;
	static double ZMAX  	= -42.;
	static double Z0      = -42.;
	static double NONLINEAR = -42;
	static int	PSNLTYPE	=-42;
	static int	EH	=-42;
	
	static double *table_p, *table_m, *thet;
	static double dtheta = .0, logthetamin = 0., logthetamax=0.,lnlc=0.,dlnl=0.;
	static int f_table=0,nc=0;
	double res, thetalog, t;
	FILE *F;
	char *xi_name;


	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN != bg_zdistrpar->zmin || ZMAX != bg_zdistrpar->zmax|| NONLINEAR != cosmopar->nonlinear || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		ZMIN = bg_zdistrpar->zmin;
		ZMAX = bg_zdistrpar->zmax;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;

		nc=N_thetaH/2+1;
		lnlc=0.5*(log(s_max)+log(s_min));
		dlnl=(log(s_max)-log(s_min))/(1.0*N_thetaH-1);
		logthetamin=(nc-(N_thetaH-1))*dlnl-lnlc;
		logthetamax=nc*dlnl-lnlc;
		dtheta=(logthetamax-logthetamin)/(1.0*N_thetaH-1.0);

		if(!f_table)
		{
			table_p=calloc(N_thetaH,sizeof(double));
			table_m=calloc(N_thetaH,sizeof(double));
			thet=calloc(N_thetaH,sizeof(double));
			f_table=1;
		}

		xi_name = (char*)malloc(200*sizeof(char));
		mkdir(path,0777);

		sprintf(xi_name, "%s/xi-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f_%d",path,cosmopar->omm,cosmopar->omv,cosmopar->Gamma,cosmopar->sigma8,cosmopar->n,bg_zdistrpar->beta,bg_zdistrpar->z0,bg_zdistrpar->zmin,bg_zdistrpar->zmax,cosmopar->nonlinear);
  //printf("Opening %s\n",xi_name);

  		//read xi from file
		if ((F=fopen(xi_name, "r")))
		{
			fscanf(F, "0.0 %lf\n", &xi_plus_zero);
			thetalog = logthetamin;
			i = 0;
			while (fscanf(F, "%lf %lf %lf\n", &t, table_p+i, table_m+i) != EOF)
			{
				if (fabs(t - exp(thetalog)) > epsilon3)
				{
					printf("xi: t:%g, tlog: %g\n",t,exp(thetalog));
					printf("Inconsistency with xi file %s!\n",xi_name); exit(1);
				}
				thetalog += dtheta;
				i++;
			}
			fileclose(F);
		}
		else
		{
   //calculate xi and write to disk
			F = fileopen(xi_name, "w");

			Hankel_forXi(1,table_p,thet);
			Hankel_forXi(0,table_m,thet);
   //xi_plus_zero = int_over_p2_j0(0.0); //should be replaced by something better...
   
			fprintf(F, "0.0 %5.20e\n", xi_plus_zero);
			thetalog = logthetamin;

   //Needs to be replaced ...
			for (i=0; i<N_thetaH; i++)
			{
				fprintf(F, "%5.10e %5.20e %5.20e\n", thet[i], table_p[i], table_m[i]);
			}
			fileclose(F);
			free((void*)xi_name);
		}
	}

 //interpolate xi
	if (theta>=theta_min)
	{
		if (theta>theta_max)
		{
			printf("theta=%e larger than theta_max=%e in xi\n", theta, theta_max);
			error("theta too large in xi");
		}
		thetalog = log(theta);
		t = (thetalog-logthetamin)/dtheta;
		i = (int)(floor(t));
		if (pm == +1)
		{
			res = (t-i)*(table_p[i+1] - table_p[i]) + table_p[i];
		}
		else
		{
			res = (t-i)*(table_m[i+1] - table_m[i]) + table_m[i];
		}
  
		return res;
	}
	else if (theta < epsilon2)
	{
		if (pm == +1) { return xi_plus_zero; }
		else { return 0.0; }
	}
	else
	{
		t = theta*exp(-logthetamin);
		if (pm == +1) { return t*(table_p[0]-xi_plus_zero) + xi_plus_zero; }
		else { return t*(table_m[0]); }
	}
}


#define hankel_r_min 0.01   // [Mpc/h] min. separation for which routine works
#define hankel_r_max 300.0  // [Mpc/h] max. separation for which routine works
#define hankel_k_min 1.e-5  // [h/Mpc] min. wavenumber 
#define hankel_k_max 1.e5   // [h/Mpc] max. wavenumber 
double Hankel_Xi_3D(int multipole,double r,double a)
{
  int i;
  double res, rlog, t;

  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0      = -42.;
  static double WA      = -42.;
  static double SIGMA_8 = -42.;
  static double GAMMA   = -42.;
  static double N_SPEC  = -42.;
  static double SCALE_A = -42.;
  static int EH	        = -42;
  static int MULTIPOLE  = -42;
  static int NONLINEAR  = -42;
	
  static double *table_xi, *rbins;
  static double dr = .0, logrmin = 0., logrmax=0.,lnlc=0.,dlnl=0.;
  static int f_table=0,nc=0;
  double prefac=sqrt(PI/2.)/(2.*PI*PI);
  if (multipole%2==0) prefac*=pow(-1.,multipole/2);  //ignore prefactors of i for odd multipoles

  if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || EH != cosmopar->transfer_EH || MULTIPOLE != multipole || NONLINEAR != cosmopar->nonlinear || SCALE_A != a) {
    OMEGA_M = cosmopar->omm;
    OMEGA_V = cosmopar->omv;
    W0      = cosmopar->w0;
    WA      = cosmopar->wa;
    SIGMA_8 = cosmopar->sigma8;
    N_SPEC  = cosmopar->n;
    GAMMA   = cosmopar->Gamma;
    EH	=	cosmopar->transfer_EH;
    MULTIPOLE = multipole;
    NONLINEAR = cosmopar->nonlinear;
    SCALE_A = a;
    
    nc=N_thetaH/2+1;
    lnlc=0.5*(log(hankel_k_max)+log(hankel_k_min));
    dlnl=(log(hankel_k_max)-log(hankel_k_min))/(1.0*N_thetaH-1);
    logrmin=(nc-(N_thetaH-1))*dlnl-lnlc;
    logrmax=nc*dlnl-lnlc;
    dr=(logrmax-logrmin)/(1.0*N_thetaH-1.0);
    
    if(!f_table) {
      table_xi=calloc(N_thetaH,sizeof(double));
      rbins=calloc(N_thetaH,sizeof(double));
      f_table=1;
    }

    Hankel_forXi_3D(multipole,a,table_xi,rbins);
  }

  //interpolate xi
  if ((r>hankel_r_max)||(r<hankel_r_min)) {
    printf("Error in routine Hankel_Xi_3D: input separation %g outside valid range [%g;%g].\n",r,hankel_r_min,hankel_r_max);
    exit(-1);
  }

  rlog = log(r);
  t = (rlog-logrmin)/dr;
  i = (int)(floor(t));
  res = (t-i)*(table_xi[i+1] - table_xi[i]) + table_xi[i];
  return(prefac*res);
}


void Hankel_forXi_3D(int multipole,double a,double *xi,double *r)
{
  int i,nc;
  double dlnk,logkmax,logkmin,k,kk,arg[2],dlnr,lnrmin,lnrmax,lnrc;
  double *k2P;
  fftw_plan plan1,plan;
  fftw_complex *f_k2P,*conv;
  fftw_complex kernel;
  const double coverh=ckms/100.;

  k2P=fftw_malloc(N_thetaH*sizeof(double));
  f_k2P=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));
  conv=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));

  plan=fftw_plan_dft_r2c_1d(N_thetaH,k2P,f_k2P,FFTW_ESTIMATE);
  plan1=fftw_plan_dft_c2r_1d(N_thetaH,conv,k2P,FFTW_ESTIMATE);

  logkmax=log(hankel_k_max);
  logkmin=log(hankel_k_min);
  dlnk=(logkmax-logkmin)/(1.0*N_thetaH-1);
  lnrmin=log(hankel_r_min); 
  lnrmax=log(hankel_r_max);
  dlnr=(lnrmax-lnrmin)/(1.0*N_thetaH-1);
  lnrc=0.5*(logkmax+logkmin);
  nc=N_thetaH/2+1;

  // power spectrum on logarithmic bins
  for(i=0;i<N_thetaH;i++) {
    k=exp(lnrc+(i-nc)*dlnk);
    if (cosmopar->nonlinear) k2P[i]=k*k*P_NL(a,k*coverh)*ps_ch_rescale;
    else k2P[i]=k*k*P_L(a,k*coverh)*ps_ch_rescale;   //sample k^2*P_delta on logarithmic bins
  }
	
  // go to log-Fourier-space
  fftw_execute(plan);
	
  arg[0]=-0.5;   // bias q
  arg[1]=2.*multipole+0.5;  // order of Bessel function

  // perform the convolution, negative sign for kernel (cc!)
  for(i=0;i<N_thetaH/2+1;i++) {
    kk=2*PI*i/(dlnk*N_thetaH);
    Hankel_Kernel_FT(kk,&kernel,arg,2);
    conv[i][0]=f_k2P[i][0]*kernel[0]-f_k2P[i][1]*kernel[1];
    conv[i][1]=f_k2P[i][1]*kernel[0]+f_k2P[i][0]*kernel[1];
  }
	
  conv[0][1]=0;	//force Nyquist- and 0-frequency-components to be double
  conv[N_thetaH/2][1]=0;

  // go back to real space, i labels log-bins in theta
  fftw_execute(plan1);
	
  for(i=0;i<N_thetaH;i++) {
    r[N_thetaH-i-1]=exp((nc-i)*dlnk-lnrc);          //r=1/k
    xi[N_thetaH-i-1]=k2P[i]/(r[N_thetaH-i-1]*N_thetaH);
  }
 
  // clean up
  fftw_free(conv);
  fftw_free(k2P);
  fftw_free(f_k2P);
  fftw_destroy_plan(plan);
  fftw_destroy_plan(plan1);
}
#undef hankel_r_min
#undef hankel_r_max
#undef hankel_k_min
#undef hankel_k_max




double int_for_wgg_rsd(double chi, void *args)
{
  double *arg=(double*)args;  //[0]: z; [1]: rp; [2]: multipole
  int multipole=(int)arg[2];
  double a=1./(1.+arg[0]);
  double r=sqrt(chi*chi+arg[1]*arg[1]);
  double mu=chi/r;
  double res=Hankel_Xi_3D(multipole,r,a)*gsl_sf_legendre_Pl(multipole,mu);
  return(res);
}


void wgg_rsd(double *rp,int nr,double z,double pi_min,double pi_max,double **res)
{
  int i,j;
  double a,f,arg[3],xil[3];

  arg[0]=z;
  a=1./(1.+z);
  f=growth_rate(a);

  gsl_set_error_handler_off();
  for (i=0;i<nr;i++) {
    arg[1]=rp[i];
    for (j=0;j<3;j++) {
      arg[2]=2*j+0.1;
      xil[j]=int_GSL_integrate_qag(int_for_wgg_rsd,arg,pi_min,pi_max,NULL,2048);
    }
    res[0][i]=xil[0];
    res[1][i]=f*(2./3.*xil[0]+4./3.*xil[1]);
    res[2][i]=f*f*(1./5.*xil[0]+4./7.*xil[1]+8./35.*xil[2]);
  }
  return;
}



void Hankel_forXi(int pm, double *xi, double *thet)
{
	double               dlnl,loglmax,loglmin,l,kk,arg[2],dlntheta,lntmin,lntmax,lnrc;
	double               *lP;
	fftw_plan            plan1,plan;
	fftw_complex 	      *f_lP,*conv;
	fftw_complex         kernel;
	int                  i,nc;

// 	char fname[100];
// 	FILE      *fp;

	lP=fftw_malloc(N_thetaH*sizeof(double));
	f_lP=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));
	conv=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));

	plan=fftw_plan_dft_r2c_1d(N_thetaH,lP,f_lP,FFTW_ESTIMATE);
	plan1=fftw_plan_dft_c2r_1d(N_thetaH,conv,lP,FFTW_ESTIMATE);

 //declare these static!
	loglmax=log(s_max);
	loglmin=log(s_min);
	dlnl=(loglmax-loglmin)/(1.0*N_thetaH-1);
	lntmin=log(theta_min); lntmax=log(theta_max);
	dlntheta=(lntmax-lntmin)/(1.0*N_thetaH-1);
	lnrc=0.5*(loglmax+loglmin);
	nc=N_thetaH/2+1;

 //Power spectrum on logarithmic bins
	for(i=0;i<N_thetaH;i++)
	{
		l=exp(lnrc+(i-nc)*dlnl);
		lP[i]=l*P_2(l);   //sample l*P_2 on logarithmic bins
	}
	

 //Go to log-Fourier-space
	
	fftw_execute(plan);
	
	
	arg[0]=0;
	arg[1]=(pm ? 0 : 4);

 //perform the convolution, negative sign for kernel (cc!)
	
	for(i=0;i<N_thetaH/2+1;i++)
	{
		kk=2*PI*i/(dlnl*N_thetaH);
		Hankel_Kernel_FT(kk,&kernel,arg,2);
		conv[i][0]=f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
		conv[i][1]=f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
	}
	
	conv[0][1]=0;			//force Nyquist- and 0-frequency-components to be double
	conv[N_thetaH/2][1]=0;

 //go back to double space, i labels log-bins in theta
	
	fftw_execute(plan1);
	
 //sprintf(fname,"hankel_%d",pm);
 //fp=openf(fname,"w");
	for(i=0;i<N_thetaH;i++)
	{
		thet[N_thetaH-i-1]=exp((nc-i)*dlnl-lnrc);          //theta=1/l
		xi[N_thetaH-i-1]=lP[i]/(thet[N_thetaH-i-1]*2*PI*N_thetaH);   //sample l*P_2 on logarithmic bins
  //fprintf(fp,"%g %g\n",thet[i],xi[i]);
	}
 
 //fclose(fp);

 //Clean up
	fftw_free(conv);
	fftw_free(lP);
	fftw_free(f_lP);
	fftw_destroy_plan(plan);
	fftw_destroy_plan(plan1);
}



//Convolution kernel for Hankel-Transform with bias q, Bessel_mu
#define LN2 0.69314718
void Hankel_Kernel_FT(double x, fftw_complex *res, double *arg, int argc)
{
 //arg[0]: q, arg[1]:mu
  fftw_complex       a1,a2,g1,g2;
  double             mu,mod,xln2,si,co,d1,d2,pref,q;

  q=arg[0];
  mu=arg[1];

  //arguments for complex cosmopar->Gamma
  a1[0]=0.5*(1.0+mu+q);
  a2[0]=0.5*(1.0+mu-q);
  a1[1]=0.5*x; 
  a2[1]=-a1[1];

  math_cdgamma(a1,&g1); 
  math_cdgamma(a2,&g2);

  xln2=x*LN2;
  si=sin(xln2); 
  co=cos(xln2);
  d1=g1[0]*g2[0]+g1[1]*g2[1]; //Re
  d2=g1[1]*g2[0]-g1[0]*g2[1]; //Im
  mod=g2[0]*g2[0]+g2[1]*g2[1];
  pref=exp(LN2*q)/mod;
  
  (*res)[0]=pref*(co*d1-si*d2);
  (*res)[1]=pref*(si*d1+co*d2);
}
#undef LN2


double Hankel_Xix(int pm, double theta, int bin1, int bin2,char *path)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double ZMIN  = -42.;
	static double ZMAX      = -42.;
	static double NONLINEAR = -42;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double *table_p, *table_m, *thet;
	static double dtheta = .0, logthetamin = 0., logthetamax = 0., lnlc=0., dlnl=0.;
	static int f_table=0,nc=0;
	double arg[2];
	double res, thetalog, t;
	int i;
	FILE *F;
	char *xi_name;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN != bg_zdistrpar->zmin || ZMAX != bg_zdistrpar->zmax|| NONLINEAR != cosmopar->nonlinear|| NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		ZMIN = bg_zdistrpar->zmin;
		ZMAX = bg_zdistrpar->zmax;
		NBIN1 = bin1;
		NBIN2 = bin2;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
		
		logthetamin = log(theta_min);
		logthetamax = log(theta_max);
		dtheta = (logthetamax - logthetamin)/(N_thetaH - 1.0);

		arg[0]=bin1;
		arg[1]=bin2;

		nc=N_thetaH/2+1;
		lnlc=0.5*(log(s_max)+log(s_min));
		dlnl=(log(s_max)-log(s_min))/(1.0*N_thetaH-1);
		logthetamin=(nc-(N_thetaH-1))*dlnl-lnlc;
		logthetamax=nc*dlnl-lnlc;
		dtheta=(logthetamax-logthetamin)/(1.0*N_thetaH-1.0);

		if(!f_table)
		{
			table_p=calloc(N_thetaH,sizeof(double));
			table_m=calloc(N_thetaH,sizeof(double));
			thet=calloc(N_thetaH,sizeof(double));
			f_table=1;
		}

      //printf("Func Xi: dtheta=%f\n",dtheta);
		xi_name = (char*)malloc(300*sizeof(char));
		mkdir(path,0777);
		sprintf(xi_name, "%s/xi_%d-%d_-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f-%.6f", path, bin1, bin2, cosmopar->omm, cosmopar->omv, cosmopar->Gamma, cosmopar->sigma8, cosmopar->n, bg_zdistrpar->beta, bg_zdistrpar->z0);
		// read xi from file 
		if ((!bg_zdistrpar->zfile)&&(F=fopen(xi_name, "r")))
		{
			fscanf(F, "0.0 %lf\n", &xi_plus_zero);
			thetalog = logthetamin;
			i = 0;
			while (fscanf(F, "%lf %lf %lf\n", &t, table_p+i, table_m+i) != EOF)
			{
				if (fabs(t - exp(thetalog)) > epsilon3)
				{
					printf("xi: t:%g, tlog: %g\n",t,exp(thetalog));
					printf("Inconsistency with xi file %s!\n",xi_name); exit(1);
				}
				thetalog += dtheta;
				i++;
			}
			fileclose(F);
		}
		else
		{
		  // calculate xi and write to disk 
			F = fileopen(xi_name, "w");
			thetalog = logthetamin;
	 
			Hankel_forXix(1,table_p,thet,bin1,bin2);
			Hankel_forXix(0,table_m,thet,bin1,bin2);

			fprintf(F, "0.0 %5.20e\n", xi_plus_zero);
			thetalog = logthetamin;
			for (i=0; i<N_thetaH; i++)
			{
				fprintf(F, "%5.10e %5.20e %5.20e\n", thet[i], table_p[i], table_m[i]);
			}
			fileclose(F);
			free((void*)xi_name);
		}
	}
// interpolate xi 
	if (theta>=theta_min)
	{
		if (theta>theta_max)
		{
			printf("theta=%e larger than theta_max=%e in xi\n", theta, theta_max);
			error("theta too large in xi");
		}
		thetalog = log(theta);
		t = (thetalog-logthetamin)/dtheta;
		i = (int)(floor(t));
		if (pm == +1)
		{
			res = (t-i)*(table_p[i+1] - table_p[i]) + table_p[i];
		}
		else
		{
			res = (t-i)*(table_m[i+1] - table_m[i]) + table_m[i];
		}
		return res;
	}
	else if (theta < epsilon2)
	{
		if (pm == +1) { return xi_plus_zero; }
		else { return 0.0; }
	}
	else
	{
		t = theta/theta_min;
		if (pm == +1) { return t*(table_p[0]-xi_plus_zero) + xi_plus_zero; }
		else { return t*(table_m[0]); }
	}
}



double Hankel_Xi_gicorr(int pm, double theta, int bin1, int bin2)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double ZMIN  = -42.;
	static double ZMAX      = -42.;
	static double NONLINEAR = -42;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static int PM=-42;
	static double *table, *thet;
	static double dtheta = .0, logthetamin = 0., logthetamax = 0., lnlc=0., dlnl=0.;
	static int f_table=0,nc=0;
	double arg[2];
	double res, thetalog, t;
	int i;
	if (pm==10) pm=1;  //recast for GG(+)
	if (pm==11) pm=0;  //recast for GG(-)

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN != bg_zdistrpar->zmin || ZMAX != bg_zdistrpar->zmax|| NONLINEAR != cosmopar->nonlinear|| NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH || PM != pm)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		ZMIN = bg_zdistrpar->zmin;
		ZMAX = bg_zdistrpar->zmax;
		NBIN1 = bin1;
		NBIN2 = bin2;
		PSNLTYPE	= cosmopar->psnltype;
		EH	= cosmopar->transfer_EH;
		PM      = pm;
		
		logthetamin = log(theta_min);
		logthetamax = log(theta_max);
		dtheta = (logthetamax - logthetamin)/(N_thetaH - 1.0);

		arg[0]=bin1;
		arg[1]=bin2;

		nc=N_thetaH/2+1;
		lnlc=0.5*(log(s_max)+log(s_min));
		dlnl=(log(s_max)-log(s_min))/(1.0*N_thetaH-1);
		logthetamin=(nc-(N_thetaH-1))*dlnl-lnlc;
		logthetamax=nc*dlnl-lnlc;
		dtheta=(logthetamax-logthetamin)/(1.0*N_thetaH-1.0);

		if(!f_table)
		{
			table=calloc(N_thetaH,sizeof(double));
			thet=calloc(N_thetaH,sizeof(double));
			f_table=1;
		}
		Hankel_forXix(pm,table,thet,bin1,bin2); 
	}

	/* interpolate xi */
	if (theta>=theta_min) {
	  if (theta>theta_max) {
	    printf("theta=%e larger than theta_max=%e in Hankel_Xi_gicorr!\n", theta, theta_max);
	    exit(-1);
	  }
	  thetalog = log(theta);
	  t = (thetalog-logthetamin)/dtheta;
	  i = (int)(floor(t));
	  res = (t-i)*(table[i+1] - table[i]) + table[i];
	  return res;
	}
	else {
	  printf("theta=%e smaller than theta_min=%e in Hankel_Xi_gicorr!\n", theta, theta_min);
	  exit(-1);
	}
}

double Hankel_Xi_gicorr_spec(int pm, double R, double z)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double ZMIN    = -42.;
	static double ZMAX    = -42.;
	static double Z       = -42.;
	static double NONLINEAR = -42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static int PM=-42;
	static double logkmin=-8.,logkmax=5.;  //determines maximum R range
	static double *table, *r;
	static double dr = .0, logrmin = 0., logrmax = 0., lnlc=0., dlnl=0.;
	static int f_table=0,nc=0;
	double res, rlog, t;
	int i;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN != bg_zdistrpar->zmin || ZMAX != bg_zdistrpar->zmax|| NONLINEAR != cosmopar->nonlinear|| PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH || PM != pm || Z != z)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		ZMIN = bg_zdistrpar->zmin;
		ZMAX = bg_zdistrpar->zmax;
		PSNLTYPE	= cosmopar->psnltype;
		EH	= cosmopar->transfer_EH;
		PM      = pm;
		Z       = z;
		
		nc=N_thetaH/2+1;
		lnlc=0.5*(logkmax+logkmin);
		dlnl=(logkmax-logkmin)/(1.0*N_thetaH-1);
		logrmin=(nc-(N_thetaH-1))*dlnl-lnlc;
		logrmax=nc*dlnl-lnlc;
		dr=(logrmax-logrmin)/(1.0*N_thetaH-1.0);

		if(!f_table)
		{
			table=calloc(N_thetaH,sizeof(double));
			r=calloc(N_thetaH,sizeof(double));
			f_table=1;
		}
		Hankel_forXix_spec(pm,table,r,z,logkmin,logkmax); 
	}

	/* interpolate w */
	if (R>=exp(logrmin)) {
	  if (R>exp(logrmax)) {
	    printf("R=%e larger than R_max=%e in Hankel_Xi_gicorr_spec!\n",R,exp(logrmax));
	    exit(-1);
	  }
	  rlog = log(R);
	  t = (rlog-logrmin)/dr;
	  i = (int)(floor(t));
	  res = (t-i)*(table[i+1] - table[i]) + table[i];
	  return res;
	}
	else {
	  printf("R=%e smaller than R_min=%e in Hankel_Xi_gicorr_spec!\n",R,exp(logrmin));
	  exit(-1);
	}
}



double Hankel_Xix_mc(int pm, double theta, int bin1, int bin2)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double SIGMA_8 = -42.;
	static double GAMMA   = -42.;
	static double N_SPEC  = -42.;
	static double BETA_P  = -42.;
	static double Z0      = -42.;
	static double ZMIN  	= -42.;
	static double ZMAX      = -42.;
	static double NONLINEAR = -42;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static double *table_p, *table_m, *thet;
	static double dtheta = .0, logthetamin = 0., logthetamax = 0., lnlc=0., dlnl=0.;
	static int f_table=0,nc=0;
	double arg[2];
	double res, thetalog, t;
	int i;

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || SIGMA_8 != cosmopar->sigma8 || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || ZMIN != bg_zdistrpar->zmax || ZMAX != bg_zdistrpar->zmin|| NONLINEAR != cosmopar->nonlinear|| NBIN1 != bin1 || NBIN2 != bin2)
	{
		OMEGA_M = cosmopar->omm;
		OMEGA_V = cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		SIGMA_8 = cosmopar->sigma8;
		N_SPEC  = cosmopar->n;
		GAMMA   = cosmopar->Gamma;
		BETA_P  = bg_zdistrpar->beta;
		Z0      = bg_zdistrpar->z0;
		NONLINEAR = cosmopar->nonlinear;
		ZMIN = bg_zdistrpar->zmin;
		ZMAX = bg_zdistrpar->zmax;
		NBIN1 = bin1;
		NBIN2 = bin2;
		logthetamin = log(theta_min);
		logthetamax = log(theta_max);
		dtheta = (logthetamax - logthetamin)/(N_thetaH - 1.0);

		arg[0]=bin1;
		arg[1]=bin2;

		nc=N_thetaH/2+1;
		lnlc=0.5*(log(s_max)+log(s_min));
		dlnl=(log(s_max)-log(s_min))/(1.0*N_thetaH-1);
		logthetamin=(nc-(N_thetaH-1))*dlnl-lnlc;
		logthetamax=nc*dlnl-lnlc;
		dtheta=(logthetamax-logthetamin)/(1.0*N_thetaH-1.0);

		if(!f_table)
		{
			table_p=calloc(N_thetaH,sizeof(double));
			table_m=calloc(N_thetaH,sizeof(double));
			thet=calloc(N_thetaH,sizeof(double));
			f_table=1;
		}

      //printf("Func Xi: dtheta=%f\n",dtheta);
      
		thetalog = logthetamin;
      
		Hankel_forXix(1,table_p,thet,bin1,bin2);
		Hankel_forXix(0,table_m,thet,bin1,bin2);
	}
	/* interpolate xi */
	if (theta>=theta_min)
	{
		if (theta>theta_max)
		{
			printf("theta=%e larger than theta_max=%e in xi\n", theta, theta_max);
			error("theta too large in xi");
		}
		thetalog = log(theta);
		t = (thetalog-logthetamin)/dtheta;
		i = (int)(floor(t));
		if (pm == +1)
		{
			res = (t-i)*(table_p[i+1] - table_p[i]) + table_p[i];
		}
		else
		{
			res = (t-i)*(table_m[i+1] - table_m[i]) + table_m[i];
		}
		return res;
	}
	else if (theta < epsilon2)
	{
		if (pm == +1) { return xi_plus_zero; }
		else { return 0.0; }
	}
	else
	{
		t = theta/theta_min;
		if (pm == +1) { return t*(table_p[0]-xi_plus_zero) + xi_plus_zero; }
		else { return t*(table_m[0]); }
	}
}



void Hankel_forXix(int pm, double *xi, double *thet, int bin1, int bin2)
{
	double               dlnl,loglmax,loglmin,l,kk,arg[2],dlntheta,lntmin,lntmax,lnrc;
	double               *lP;
	fftw_plan            plan1,plan;
	fftw_complex 	      *f_lP,*conv;
	fftw_complex         kernel;
	int                  i,nc;

	lP=fftw_malloc(N_thetaH*sizeof(double));
	f_lP=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));
	conv=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));

	plan=fftw_plan_dft_r2c_1d(N_thetaH,lP,f_lP,FFTW_ESTIMATE);
	plan1=fftw_plan_dft_c2r_1d(N_thetaH,conv,lP,FFTW_ESTIMATE);

 //declare these static!
	loglmax=log(s_max);
	loglmin=log(s_min);
	dlnl=(loglmax-loglmin)/(1.0*N_thetaH-1);
	lntmin=log(theta_min); lntmax=log(theta_max);
	dlntheta=(lntmax-lntmin)/(1.0*N_thetaH-1);
	lnrc=0.5*(loglmax+loglmin);
	nc=N_thetaH/2+1;

 //Power spectrum on logarithmic bins
	if (pm<=1) {        //GG
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_2x(l,bin1,bin2);   //sample l*P_2 on logarithmic bins
	  }
	}
	else if (pm==2) {   //gg
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_2pp(l,bin1,bin2,1);
	  }
	}
	else if (pm==3) {   //gG
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_2pq(l,bin1,bin2,1);
	  }
	}
	else if (pm==4) {   //gI
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_2pp(l,bin1,bin2,0);
	  }
	}
	else if (pm==5) {   //mG
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_2x(l,bin1,bin2);
	  }
	}
	else if ((pm==6)||(pm==7)) {   //II
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_2IIx(l,bin1,bin2);
	  }
	}
	else if ((pm==8)||(pm==9)) {   //GI
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_2iax(l,bin1,bin2);
	  }
	}
	else {
	  printf("Error: wrong argument pm in Hankel_forXix!\n");
	  exit(-1);
	} 

 //Go to log-Fourier-space
	fftw_execute(plan);

	arg[0]=0;
	//arg[1]=(pm ? 0 : 4);  //old version
	if (pm==0) arg[1]=4;   //xi_-
	if (pm==1) arg[1]=0;   //xi_+
	if (pm==2) arg[1]=0;   //gg  
	if (pm==3) arg[1]=2;   //gG  
	if (pm==4) arg[1]=2;   //gI  
	if (pm==5) arg[1]=2;   //mG 
	if (pm==6) arg[1]=0;   //II(+)
	if (pm==7) arg[1]=4;   //II(-)
	if (pm==8) arg[1]=0;   //GI(+)
	if (pm==9) arg[1]=4;   //GI(-)

 //perform the convolution, negative sign for kernel (cc!)
	for(i=0;i<N_thetaH/2+1;i++)
	{
		kk=2*PI*i/(dlnl*N_thetaH);
		Hankel_Kernel_FT(kk,&kernel,arg,2);
		conv[i][0]=f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
		conv[i][1]=f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
	}
	conv[0][1]=0;			//Nyquist- and 0-frequency-components are double!
	conv[N_thetaH/2][1]=0;

 //go back to double space, i labels log-bins in theta
	fftw_execute(plan1);

	for(i=0;i<N_thetaH;i++)
	{
		thet[N_thetaH-i-1]=exp((nc-i)*dlnl-lnrc);
		xi[N_thetaH-i-1]=lP[i]/(thet[N_thetaH-i-1]*2*PI*N_thetaH) ;   //sample l*P_2 on logarithmic bins
	}

 //Clean up
	fftw_free(conv);
	fftw_free(lP);
	fftw_free(f_lP);
	fftw_destroy_plan(plan);
	fftw_destroy_plan(plan1);
}


void Hankel_forXix_spec(int pm, double *xi, double *r, double z, double logkmin, double logkmax)
{
	double               dlnl,l,kk,arg[2],dlnr,lnrc,a;
	double               *lP;
	fftw_plan            plan1,plan;
	fftw_complex 	     *f_lP,*conv;
	fftw_complex         kernel;
	int                  i,nc;
	const double coverh=ckms/100.;

	lP=fftw_malloc(N_thetaH*sizeof(double));
	f_lP=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));
	conv=fftw_malloc((N_thetaH/2+1)*sizeof(fftw_complex));

	plan=fftw_plan_dft_r2c_1d(N_thetaH,lP,f_lP,FFTW_ESTIMATE);
	plan1=fftw_plan_dft_c2r_1d(N_thetaH,conv,lP,FFTW_ESTIMATE);

	//declare these static!
	dlnl=(logkmax-logkmin)/(1.0*N_thetaH-1);
	//dlnr=(logrmax-logrmin)/(1.0*N_thetaH-1);
	lnrc=0.5*(logkmax+logkmin);
	nc=N_thetaH/2+1;

	//Power spectrum on logarithmic bins
	a=1./(1.+z);
	if (pm==0) {   //gg
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_NL(a,l*coverh)*ps_ch_rescale;  
	  }
	}
	else if (pm==2) {   //gI
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*(-1.)*P_GI(a,l*coverh)*ps_ch_rescale; //turn it positive
	  }
	}
	else if ((pm==4)||(pm==5)) {   //II
	  for(i=0;i<N_thetaH;i++) {
	    l=exp(lnrc+(i-nc)*dlnl);
	    lP[i]=l*P_II(a,l*coverh)*ps_ch_rescale; 
	  }
	}
	else {
	  printf("Error: wrong argument pm in Hankel_forXix_spec!\n");
	  exit(-1);
	} 

	//Go to log-Fourier-space
	fftw_execute(plan);

	arg[0]=0;
	if (pm==0) arg[1]=0;
	if (pm==2) arg[1]=2;
	if (pm==4) arg[1]=0;
	if (pm==5) arg[1]=4;
	if (pm==8) arg[1]=0;
	if (pm==9) arg[1]=4;

	//perform the convolution, negative sign for kernel (cc!)
	for(i=0;i<N_thetaH/2+1;i++)
	{
		kk=2*PI*i/(dlnl*N_thetaH);
		Hankel_Kernel_FT(kk,&kernel,arg,2);
		conv[i][0]=f_lP[i][0]*kernel[0]-f_lP[i][1]*kernel[1];
		conv[i][1]=f_lP[i][1]*kernel[0]+f_lP[i][0]*kernel[1];
	}
	conv[0][1]=0;  //Nyquist- and 0-frequency-components are double!
	conv[N_thetaH/2][1]=0;

	//go back to double space, i labels log-bins in theta
	fftw_execute(plan1);

	for(i=0;i<N_thetaH;i++)
	{
	  r[N_thetaH-i-1]=exp((nc-i)*dlnl-lnrc);
	  xi[N_thetaH-i-1]=lP[i]/(r[N_thetaH-i-1]*2*PI*N_thetaH) ;   //sample l*P_2 on logarithmic bins
	}

	//Clean up
	fftw_free(conv);
	fftw_free(lP);
	fftw_free(f_lP);
	fftw_destroy_plan(plan);
	fftw_destroy_plan(plan1);
}





// function for variable dark energy equation of state
double omv_vareos(double a)
{
  return(cosmopar->omv*exp(-3.*((cosmopar->w0+cosmopar->wa+1)*log(a)+cosmopar->wa*(1-a))));
}




// 3D GI power spectrum, 1-halo term fit                                                                                                                 
double P_GI_1h(double a, double k)
{
  double z=1./a-1.;
  double p1=0.02056*exp(5.909*pow(z,0.3798));  // in h/Mpc                                                                                                
  double p2=1.978*exp(1.087*pow(z,0.6655));
  double p3=4.154*exp(0.1912*pow(z,0.4368));
  return (-1.)*gammascale_sb10*pow(k/p1,2.)/(1.+pow(k/p2,p3));
}


// 3D IA power spectrum 
double P_GI(double a, double k)
{
  double ps;
  //const double C1rhocr=-0.0134; //C_1*rho_cr  -- now defined in ps.h
  double fac=C1rhocr*cosmopar->omm*IAPS_MACRO(a)/(growfac(a,1.,1.)/growfac(1.,1.,1.));
  const double coverh=ckms/100.;

  if (cosmopar->iamethod==0) ps=fac*P_L(a,k);
  else if (cosmopar->iamethod==1) ps=fac*P_NL(a,k);
  else if (cosmopar->iamethod==2) ps=fac*P_L(a,k)+P_GI_1h(a,k/coverh)*ps_ch_rescale;
  else {
    printf("Error in routine P_GI: no method to compute power spectrum.\n");
    exit(-1);
  }
  return ps;
}


// Toy model (power law) 3D IA power spectrum
double P_TOY(double a, double k)
{
  const double slope=cosmopar->toyexp;   // def 0.4 
  const double frac=1.e0;   // def 1.
  const double kref=1.e0;  // def 1.
  const double aref=1./(1.+bg_zdistrpar->z0*1.412);
  const double coverh=ckms/100.;

  double amp=frac*fabs(P_GI(aref,kref*coverh));
  return amp*pow(k/(kref*coverh),slope-2.)*pow(a,-1.);  //def -1.
}


// matter power spectrum + galaxy bias model for gm and gg signals 
double P_bias(double a, double k)
{
        const double offset=cosmopar->toyexp;   // b_g=(z+offset)
	double fac,z=1./a-1.;
	if (offset==42.) fac=1.;  // unit bias
	else {
	  if (cosmopar->ia==3) fac=(z+offset); // for mg	
	  if (cosmopar->ia==2) fac=(z+offset)*(z+offset); //for gg
	}
	return fac*P_NL(a,k);
}


// 2D GI power spectrum
double int_for_p_2iax(double a, void *args)
{
	double hoverh0, asqr, s, fKw, f, res;
	int bin1,bin2;

	double *arg	=	(double*)args;
		
	bin1=(int)arg[0];
	bin2=(int)arg[1];

	if (a >= 1.0) error("a>=1 in int_for_p_2iax");
	s       = sglob;
	asqr    = a*a;
	fKw     = f_K(w(a));
	f       = s/fKw;
	//hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));
	res = (g_sourcex1(a,bin1)*probx2(1./a-1.,bin2)+probx1(1./a-1.,bin1)*g_sourcex2(a,bin2))/(asqr*a*fKw);  // standard symmetrised

	//res = probx1(1./a-1.,bin1)*g_sourcex2(a,bin2)/(asqr*a*fKw); //one-sided IG
	//res = g_sourcex1(a,bin1)*probx2(1./a-1.,bin2)/(asqr*a*fKw); //one-sided GI

	if (cosmopar->toy==1) res *= P_TOY(a, f);
	if (cosmopar->toy==2) res *= P_bias(a, f);
	if (cosmopar->toy==0) res *= (-1)*P_GI(a, f); //need to turn it positive
	return res;
}


//tomography IA power spectra
double P_2iax(double s, int bin1, int bin2)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double N_SPEC  = -42.;
	static double GAMMA   = -42.;
	static double BETA_P  = -42.;
	static double ZMAX    = -42.;
	static double ZMIN    = -42.;
	static double Z0      = -42.;
        static double FCAT    = -42.;
	static double SIG     = -42.;
	static double DELTAZ  = -42.;
	static double NONLINEAR = -42.;
	static double SIGMA_8 = -42.;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table[N_s];
	static double ds = .0, logsmin = .0, logsmax = .0;
	double   ss, slog, f1, f2;
	int    i,j;
	double arg[2],integral,step,stepsize;
	#define NLOOP 10
	double aa;
	double da=(1.-a_min)/(1.*NLOOP);

	double a1 = 1./(1.+bg_zdistrpar->zmed[bin1]); 
	double a2 = 1./(1.+bg_zdistrpar->zmed[bin2]);
	double w1=w(a1);
	double w2=w(a2);


	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || NONLINEAR != cosmopar->nonlinear || SIGMA_8 != cosmopar->sigma8 || ZMAX != bg_zdistrpar->zmax || ZMIN != bg_zdistrpar->zmin || NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{	
		logsmin = log(s_min);
		logsmax = log(s_max);
		if(N_s==1) ds=0.; else
		ds = (logsmax - logsmin)/(N_s - 1.);
		slog = logsmin;
		arg[0]=bin1+0.1;
		arg[1]=bin2+0.1;

		for (i=0; i<N_s; i++, slog+=ds) {
			ss = exp(slog);
			sglob = ss;
			if (bg_zdistrpar->sheet) {
			  integral=f_K(w2-w1)/f_K(w2)/(a1*w1); //delta-function in z, not a! -> 1/a instead of 1/a^3
			  if (cosmopar->toy==1) integral *= P_TOY(a1,ss/f_K(w1));
			  if (cosmopar->toy==2) integral *= P_bias(a1,ss/f_K(w1));
			  if (cosmopar->toy==0) integral *= (-1)*P_GI(a1,ss/f_K(w1)); 
			}
			else {
			  integral=0.0;
			  for (j=0;j<NLOOP;j++) {
			    aa=a_min+j*da;
			    integral+=int_GSL_integrate_qag(int_for_p_2iax,arg,aa,aa+da,NULL,2048);
			  }
			  //integral=int_GSL_integrate_qag(int_for_p_2iax,arg,a_min,1.,NULL,2048);	
			}
			if (integral<=0.0) integral=1.e-100; //avoids nans!
			table[i] = log(3./2.*cosmopar->omm*integral); 
		}
		OMEGA_M =   cosmopar->omm;
		OMEGA_V =   cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		N_SPEC  =   cosmopar->n;
		GAMMA   =   cosmopar->Gamma;
		BETA_P  =   bg_zdistrpar->beta;
		Z0      =   bg_zdistrpar->z0;
		FCAT    =   bg_zdistrpar->photo_fcat;
		SIG     =   bg_zdistrpar->photo_sig;
		DELTAZ  =   bg_zdistrpar->photo_deltaz;
		NONLINEAR = cosmopar->nonlinear;
		SIGMA_8 =  cosmopar->sigma8;
		ZMAX    =  bg_zdistrpar->zmax;
		ZMIN    =  bg_zdistrpar->zmin;
		NBIN1 = bin1;
		NBIN2 = bin2;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
	}
	slog = log(s);
	f1 = interpol(table, N_s, logsmin, logsmax, ds, slog, cosmopar->n, cosmopar->n-4.0);
	return exp(f1);
}




// 3D II power spectrum, 1-halo term fit
double P_II_1h(double a, double k)
{
  double z=1./a-1.;
  double p1=0.01291*exp(6.781*pow(z,0.203));  // in h/Mpc
  double p2=1.98*exp(1.033*pow(z,0.7593)); 
  double p3=6.064*exp(0.1054*pow(z,0.5937)); 
  return gammascale_sb10*gammascale_sb10*pow(k/p1,4.)/(1.+pow(k/p2,p3));
}


// 3D II power spectrum
double P_II(double a, double k)
{
  double ps;
  //const double C1rhocr=0.0134; //C_1*rho_cr, now defined in ps.h
  double fac=C1rhocr*cosmopar->omm*IAPS_MACRO(a)/(growfac(a,1.,1.)/growfac(1.,1.,1.));
  const double coverh=ckms/100.;

  if (cosmopar->iamethod==0) ps=fac*fac*P_L(a,k);
  else if (cosmopar->iamethod==1) ps=fac*fac*P_NL(a,k);
  else if (cosmopar->iamethod==2) ps=fac*fac*P_L(a,k)+P_II_1h(a,k/coverh)*ps_ch_rescale;
  else {
    printf("Error in routine P_II: no method to compute power spectrum.\n");
    exit(-1);
  }
  return ps;
}


// 2D II power spectrum
double int_for_p_2IIx(double a, void *args)
{
  double asqr, s, fKw, f, res, hoverh0;
    int bin1,bin2;

    double *arg =   (double*)args;
    bin1=(int)arg[0];
    bin2=(int)arg[1];
	    
    if (a >= 1.0) error("a>=1 in int_for_p_2IIx");
    s = sglob;
    asqr = a*a;
    fKw = f_K(w(a));
    f = s/fKw;

    hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));
    res = probx1(1./a-1.,bin1)*probx2(1./a-1.,bin2)/(fKw*fKw*asqr)*hoverh0;
    if (cosmopar->toy==0) res*=P_II(a,f);
    if (cosmopar->toy==2) res*=P_bias(a,f);
    if (cosmopar->toy==1) {
      printf("Error: toy model for II power spectra currently not implemented!\n");
      exit(-1);
    }
    return res;
}

//tomography II power spectra
double P_2IIx(double s, int bin1, int bin2)
{
	static double OMEGA_M = -42.;
	static double OMEGA_V = -42.;
	static double W0      = -42.;
	static double WA      = -42.;
	static double N_SPEC  = -42.;
	static double GAMMA   = -42.;
	static double BETA_P  = -42.;
	static double ZMAX    = -42.;
	static double ZMIN    = -42.;
	static double Z0      = -42.;
        static double FCAT    = -42.;
	static double SIG     = -42.;
	static double DELTAZ  = -42.;
	static double NONLINEAR = -42.;
	static double SIGMA_8 = -42.;
	static int NBIN1=-42;
	static int NBIN2=-42;
	static int PSNLTYPE=-42;
	static int EH=-42;
	static double table[N_s];
	static double ds = .0, logsmin = .0, logsmax = .0;
	double   ss, slog, f1, f2;
	int    i,j;
	double arg[2],integral,step,stepsize;
        #define NLOOPII 10
	double aa;
	double da=(1.-a_min)/(1.*NLOOPII);

	if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || NONLINEAR != cosmopar->nonlinear || SIGMA_8 != cosmopar->sigma8 || ZMAX != bg_zdistrpar->zmax || ZMIN != bg_zdistrpar->zmin || NBIN1 != bin1 || NBIN2 != bin2 || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH)
	{
		logsmin = log(s_min);
		logsmax = log(s_max);
		ds = (logsmax - logsmin)/(N_s - 1.);
		slog = logsmin;
		arg[0]=bin1+0.1;
		arg[1]=bin2+0.1;

		for (i=0; i<N_s; i++, slog+=ds) {
			ss = exp(slog);
			sglob = ss;
			integral=0.0;
			for (j=0;j<NLOOPII;j++) {
			  aa=a_min+j*da;
			  integral+=int_GSL_integrate_qag(int_for_p_2IIx,arg,aa,aa+da,NULL,2048);
			}
  //integral=int_GSL_integrate_qag(int_for_p_2IIx,arg,a_min,1.,NULL,2048);
			table[i] = log(integral); 
		}
		OMEGA_M =   cosmopar->omm;
		OMEGA_V =   cosmopar->omv;
		W0      = cosmopar->w0;
		WA      = cosmopar->wa;
		N_SPEC  =   cosmopar->n;
		GAMMA   =   cosmopar->Gamma;
		BETA_P  =   bg_zdistrpar->beta;
		Z0      =   bg_zdistrpar->z0;
		FCAT    =   bg_zdistrpar->photo_fcat;
		SIG     =   bg_zdistrpar->photo_sig;
		DELTAZ  =   bg_zdistrpar->photo_deltaz;
		NONLINEAR = cosmopar->nonlinear;
		SIGMA_8 =  cosmopar->sigma8;
		ZMAX    =  bg_zdistrpar->zmax;
		ZMIN    =  bg_zdistrpar->zmin;
		NBIN1 = bin1;
		NBIN2 = bin2;
		PSNLTYPE	=	cosmopar->psnltype;
		EH	=	cosmopar->transfer_EH;
	}
	slog = log(s);
	f1 = interpol(table, N_s, logsmin, logsmax, ds, slog, cosmopar->n, cosmopar->n-4.0);
	return exp(f1);
}




//compute median redshift of bins
double int_for_zmedian(double z, void *args)
{
  int *bin=(int*)args;
  return(probx1(z,*bin));
}

void zmedian(int nbin)
{
  #include <gsl/gsl_errno.h>

  int i;
  static int flag=0;
  double inttot,intmed,med,upper,lower,check;
  int par[1];

  gsl_set_error_handler_off();
  for (i=0;i<nbin;i++) {
    par[0]=i;
    check=1.;
    med=(bg_zdistrpar->tomobin[i]+bg_zdistrpar->tomobin[i+1])/2.;
    if (((bg_zdistrpar->tomobin[i+1]-bg_zdistrpar->tomobin[i])>0.05)||(bg_zdistrpar->photo_sig>0.03)) {  // for small bins: median ~ mean
      lower=bg_zdistrpar->tomobin[i]-2.*bg_zdistrpar->photo_sig*(1+bg_zdistrpar->tomobin[i]);
      if (lower<bg_zdistrpar->zmin) lower=bg_zdistrpar->zmin;
      upper=bg_zdistrpar->tomobin[i+1]+2.*bg_zdistrpar->photo_sig*(1+bg_zdistrpar->tomobin[i+1]);
      if (upper>bg_zdistrpar->zmax) upper=bg_zdistrpar->zmax;
      
      while (fabs(check)>1.e-5) {
	intmed=int_GSL_integrate_qag(int_for_zmedian,par,bg_zdistrpar->zmin,med,NULL,2048);
	check=intmed-0.5;
	//printf("%i  %g  %g  %g  %g\n",i,check,lower,med,upper);
	if (check>0) {
	  upper=med;
	  med=(lower+med)/2.;
	}
	else {
	  lower=med;
	  med=(med+upper)/2.;
	}
      }
    }
    else {
      if (!flag) printf("Warning: Use mean instead of median due to small bin size %g %g!\n",bg_zdistrpar->tomobin[i+1]-bg_zdistrpar->tomobin[i],bg_zdistrpar->photo_sig);
      flag=1;
    }
    bg_zdistrpar->zmed[i]=med;
  }
}                               


// compute mean redshift of bins
double int_for_zmean(double z, void *args)
{
  int *bin=(int*)args;
  return(probx1(z,*bin)*z);
}

void zmean(int nbin)
{
  #include <gsl/gsl_errno.h>

  int i;
  int par[1];
  gsl_set_error_handler_off();

  for (i=0;i<nbin;i++) {
    par[0]=i;
    if(!bg_zdistrpar->photo_sig) //ohne photo-z Fehler
      bg_zdistrpar->zmean[i]=int_GSL_integrate_qag(int_for_zmean,par,bg_zdistrpar->tomobin[i],bg_zdistrpar->tomobin[i+1],NULL,2048);
    else 
      bg_zdistrpar->zmean[i]=int_GSL_integrate_qag(int_for_zmean,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);
  }			 
}


// compute correction of number densities per bin
double ncorrection(int bin)
{
  double par[1]={bin+0.1};  

  if (bg_zdistrpar->photo_sig) return(int_GSL_integrate_qag(int_for_prob_photoz1,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,1024))/int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,1024);
  else return(int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->tomobin[bin],bg_zdistrpar->tomobin[bin+1],NULL,1024))/int_GSL_integrate_qag(int_for_prob,NULL,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,1024);
}


// read arbitrary redshift distributions
double int_for_probzarb(double z, void *args)
{
  double *arg=(double*)args;
  int bin=(int) arg[0];

  int i,j;
  double var;
  FILE *dat;

  static int flag=0;
  static gsl_interp_accel *acc[50];   // => N_z <= 50!
  static gsl_interp *linear[50];   

  static double zf[10000];
  static double pf[51][10000];


  if (!flag) {
    if ((dat=fopen(bg_zdistrpar->zfile,"r"))==NULL) {
      printf("Couldn't open file %s in routine 'probzarb'\n",bg_zdistrpar->zfile);
      exit(-1);
    }
    i=0;
    j=0;
    while(fscanf(dat,"%lf",&var)!=EOF) {
      if (i==0) zf[j]=var;
      else pf[i-1][j]=var;
      i++;
      if (i==nbin+1) {              //nbin is global variable
	i=0;
	j++;
      }
      if (j==10000) {
	printf("Error in routine 'probzarb': maximum number of lines in file %s exceeded!\n",bg_zdistrpar->zfile);
	exit(-1);
      }
    }
    fclose(dat);
    if (i) {
      printf("Error in routine 'probzarb': Incorrect number of entries in file %s!\n",bg_zdistrpar->zfile);
      exit(-1);
    }
    const int dim=j;
    if (bg_zdistrpar->zmin!=zf[0]) {
      printf("Error in routine 'probzarb': Lowest redshift in %s is not zmin=%g!\n",bg_zdistrpar->zfile,bg_zdistrpar->zmin);
      exit(-1);
    }
    if (bg_zdistrpar->zmax!=zf[dim-1]) {
      printf("Error in routine 'probzarb': Highest redshift in %s is not zmax=%g!\n",bg_zdistrpar->zfile,bg_zdistrpar->zmax);
      exit(-1);
    }
 
    for (i=0;i<nbin;i++) {              //interpolation
      linear[i]=gsl_interp_alloc(gsl_interp_linear,dim);
      acc[i]=gsl_interp_accel_alloc();
      gsl_interp_init(linear[i],zf,pf[i],dim);
    }

    flag=1;
  }

  return(gsl_interp_eval(linear[bin],zf,pf[bin],z,acc[bin]));
}


double probzarb(double z,int bin)
{
  #include <gsl/gsl_errno.h>

  int i;
  double par[1];
  static double norm[50];
  static int flag=0;

  if (!flag) {               //normalization
    gsl_set_error_handler_off();

    for (i=0;i<nbin;i++) {
      par[0]=i+0.1;
      norm[i]=int_GSL_integrate_qag(int_for_probzarb,par,bg_zdistrpar->zmin,bg_zdistrpar->zmax,NULL,2048);
    }
    flag=1;
  }
  par[0]=bin+0.1;
  return(int_for_probzarb(z,par)/norm[bin]);
}


/*************************************************************************************/
/*  Functions for scoccimarro bispectrum      Takada&Jain05                                                                          */
//added by xun
/*************************************************************************************/

double B_3dbs(double a, double k1, double k2, double k3, double k_NL_bs1) 

// 3d-matter bispectrum 
//k in hMpc-1, B_3dbs in (h-1Mpc)^6
//k_NL_bs1 = k_NL_bs(a); for efficiency
{
//call F2_eff, P_NL
//call func_a, func_b, func_c: (a, k)

  if (cosmopar->bstype!=1||cosmopar->bstype!=2) 
  { 
    printf("Error in routine 'B_3dbs': currently only methods 'SC01' and 'GM12' implemented.\n");
    exit(-1);
  }

	double var, ang_k12, ang_k23, ang_k31, q1, q2, q3;
	double n_eff1, n_eff2, n_eff3;
	double F2_eff123, F2_eff231, F2_eff312;
	double func_a1, func_a2, func_a3; 
	double func_b1, func_b2, func_b3; 
	double func_c1, func_c2, func_c3; 
	double P_NL1, P_NL2, P_NL3;
	const double coverh=ckms/100.;
	double aparams[9];
	
	if(cosmopar->bstype!=1)
	{
         aparams = {0.25,3.5,0.2,1.0,2.0,-0.2,1.0,0.0,0.0};
        }
	if(cosmopar->bstype!=2)
	{
         aparams = {0.484,3.730,-0.849,0.392,1.013,-0.575,0.128,-0.722,-0.926;
        }

 
	n_eff1 = n_eff_bs(k1);
	n_eff2 = n_eff_bs(k2);
	n_eff3 = n_eff_bs(k3);

	q1 = k1/k_NL_bs1;	
	q2 = k2/k_NL_bs1;
	q3 = k3/k_NL_bs1;

	ang_k12 = -1./2.*(k1/k2+k2/k1-k3*k3/k1/k2);
	ang_k23 = -1./2.*(k2/k3+k3/k2-k1*k1/k2/k3);
	ang_k31 = -1./2.*(k3/k1+k1/k3-k2*k2/k3/k1);

	func_a1 = func_a(a, k1, n_eff1, q1);
	func_a2 = func_a(a, k2, n_eff2, q2);
	func_a3 = func_a(a, k3, n_eff3, q3);

	func_b1 = func_b(a, k1, n_eff1, q1);
	func_b2 = func_b(a, k2, n_eff2, q2);
	func_b3 = func_b(a, k3, n_eff3, q3);

	func_c1 = func_c(a, k1, n_eff1, q1);
	func_c2 = func_c(a, k2, n_eff2, q2);
	func_c3 = func_c(a, k3, n_eff3, q3);

	F2_eff123 = 5./7.*func_a1*func_a2 +1./2.*(k1/k2+k2/k1)*ang_k12*func_b1*func_b2 +2./7.*ang_k12*ang_k12*func_c1*func_c2;

	F2_eff231 = 5./7.*func_a2*func_a3+1./2.*(k2/k3+k3/k2)*ang_k23*func_b2*func_b3 +2./7.*ang_k23*ang_k23*func_c2*func_c3;

	F2_eff312 = 5./7.*func_a3*func_a1+1./2.*(k3/k1+k1/k3)*ang_k31*func_b3*func_b1 +2./7.*ang_k31*ang_k31*func_c3*func_c1;

	P_NL1 = P_NL(a, k1*coverh)*ps_ch_rescale;// P_NL1: dimension length^3, in units of Mpc^3
	P_NL2 = P_NL(a, k2*coverh)*ps_ch_rescale;
	P_NL3 = P_NL(a, k3*coverh)*ps_ch_rescale;

	var = F2_eff123*P_NL1*P_NL2 + F2_eff231*P_NL2*P_NL3 + F2_eff312*P_NL3*P_NL1;
	var = 2.*var;

	return var;
}


double func_a(double a, double k, double n_eff_BS, double qq)	//n_eff_BS is quantity, n_eff_bs is func.
{
//call k_NL_bs(a), n_eff_bs(k), growfac(a) 
//cosmopar->sigma8

	double var, sigma8_z;
	double Q3;
	double q;	//q=k/k_NL_bs
	double	omm,omv;
        double aparams[9];
	q = qq;
	omm = om_m(a);
	omv = om_v(a);

	Q3 = (4.-pow(2.,n_eff_BS))/(1.+pow(2.,(n_eff_BS+1.)));
	sigma8_z = cosmopar->sigma8*growfac(a,omm,omv)/growfac(1.0,cosmopar->omm,OMV_MACRO(1.));


	var = 1. + pow(sigma8_z,(-aparams[5]))*pow((0.7*Q3),0.5)*pow((q*aparams[2]),(n_eff_BS+aparams[3]);
	var = var/(1.+pow((q*aparams[0]),(n_eff_BS+aparams[1])));

	return var;
}

double func_b(double a, double k, double n_eff_BS, double qq)
{
//call k_NL_bs(a), n_eff_bs(k)

	double var;
	double q;	//q=k/k_NL_bs
	double aparams[9];

	q = qq;
	var = 1. + aparams[2]*(n_eff_BS+3.)*pow(q*aparams[6],(n_eff_BS+3.+aparams[7]));
	var = var/(1.+pow(q,(n_eff_BS+3.5+aparams[8])));

	return var;
}

double func_c(double a, double k, double n_eff_BS, double qq)
{
//call k_NL_bs(a), n_eff_bs(k)

	double var;
	double q;	//q=k/k_NL_bs
	double aparams[9];

	q = qq;

	var = 1. + 4.5*aparams[3]/(1.5+pow((n_eff_BS+3.),4.))*pow((aparams[4]*q),(n_eff_BS+3.));
	var = var/(1.+pow((aparams[4]*q),(n_eff_BS+3.5+aparams[8])));

	return var;
}


double f_P_L(double logk, void * params)	//ln ( the part of P_L(k*3000.) that depends on k, which is the transfer function) 
{
//k in hMpc-1

//call Tsqr(k), Delta_L_smith(k), sigma_8_sqr()
  double t, k, kk;
  const double coverh=ckms/100.;

  k = exp(logk);
  kk = k*coverh;	// kk in Hubble length-1

  if((cosmopar->transfer_EH == 1) && !cosmopar->transfer_tabulated) {
    t = Tsqr_EH(kk);
  }
  else if((cosmopar->transfer_EH == 2) && !cosmopar->transfer_tabulated) {
    t = Tsqr_EH_wiggle(kk);
  }
  else if((cosmopar->transfer_EH == 1) && cosmopar->transfer_tabulated) {
    t = Tsqr_tabulated(kk);
  }
  else if (cosmopar->transfer_EH == 0) {  //BBKS
    t = Tsqr(kk);
  }
  else if (cosmopar->transfer_EH == 3) {  //EBW
    return(Delta_L_smith(k)/(kk*kk*kk)*(2.*PI*PI));
  }
  else {
    printf("Error in routine f_P_L: no valid transfer function.\n");
    exit(-1);
  }
  return log(t);
}

double n_eff_bs(double k)	//n_eff=dln P_L/dlnk
{
//k in hMpc-1

//call double f_P_L(double k)

	gsl_function F;
	double result, abserr;
	double logk;

	logk = log(k);
     
	F.function = &f_P_L;
	F.params = 0;
     
//	printf ("f(x) = x^(3/2)\n");
     
//	printf("  Starting derivation for n_eff_bs(%f)\n", k);

	gsl_deriv_central (&F, logk, 1e-4, &result, &abserr);	//1e-5?
//	printf ("x = 2.0\n");
//	printf ("f'(x) = %.10f +/- %.10f\n", result, abserr);
//	printf ("exact = %.10f\n\n", 1.5 * sqrt(2.0));
//	printf("  Done.                                         \n");
//	printf("n_eff=%f\n", result);
	return result;
}

double k_NL_bs(double a)
{
//call double P_L(double a, double k)
//defined in ps.c:
//const double k_minbs = 0.001;
//const double k_maxbs = 1.e8;

	int i;
	int niter=0;
	double knl=0.,logknl=0.;
	double logkmin = 0., logkmax = 0., logkmid = 0.;
	double delta_L, delta_L_min, delta_L_max;
	const double coverh=ckms/100.;
	logkmin = log(k_minbs);	//k_minbs, k_maxbs static, others all variable
	logkmax = log(k_maxbs);

	delta_L_min = pow(k_minbs,3.)*P_L(a, k_minbs*coverh)*ps_ch_rescale/(2.*pi_sqr);
	delta_L_max = pow(k_maxbs,3.)*P_L(a, k_maxbs*coverh)*ps_ch_rescale/(2.*pi_sqr);

	if(delta_L_min>1. || delta_L_max<1.)
	{
	  printf("k_min & k_max not properly given: %g>1 or %g<1\n",delta_L_min,delta_L_max); 
	  exit(-99);
	}
	else	//begin searching for knl
	{
		knl = k_minbs;
//		printf("  Starting searching for k_nl(a = %f)\n", a);	
		delta_L = pow(knl,3.)*P_L(a, knl*coverh)*ps_ch_rescale/(2.*pi_sqr);

		while (fabs(delta_L-1.)>=1.e-5)	{
		  logkmid=(logkmin+logkmax)/2.;
		  if ((pow(exp(logkmid),3.)*P_L(a, exp(logkmid)*coverh)*ps_ch_rescale/(2.*pi_sqr))>1.) {
		    logkmax=logkmid;
		  }
		  else {
		    logkmin=logkmid;
		    knl=exp(logkmin);
		    delta_L = pow(knl,3.)*P_L(a, knl*coverh)*ps_ch_rescale/(2.*pi_sqr);	
		  }
		  niter+=1;
		}
	return knl;
	}
}

/*************************************************************************************/
/*  Functions for lensing tomography bispectrum                                                                                 */
/*************************************************************************************/

/* ============================================================ *
 * Global variables needed in integrand functions int_for_b_2	*					*
 * ============================================================ */
double sglobi,sglobj,sglobk;


double B_2x(double s1, double s2, double s3, int bin1, int bin2, int bin3)	//lensing tomo B_2
{
  double arg[3], integral;
  const double bispec_minell=5.e-3;

  assert(s1>=bispec_minell);
  assert(s1<=s_max);
  assert(s2>=bispec_minell);
  assert(s2<=s_max);
  assert(s3>=bispec_minell);
  assert(s3<=s_max);

  arg[0]=bin1+0.1;
  arg[1]=bin2+0.1;
  arg[2]=bin3+0.1;
  
  sglobi = s1;
  sglobj = s2;
  sglobk = s3;

  integral=int_GSL_integrate_qag(int_for_b_2x,arg,a_min,1.0,NULL,2048);
  return (27./8.*cosmopar->omm*cosmopar->omm*cosmopar->omm*integral);
}


double int_for_b_2x(double a, void *args)
{
  double hoverh0, asqr, s1,s2,s3, kai, f1,f2,f3, res, k_NL_bs1;
  int bin1,bin2,bin3;

  double *arg	=	(double*)args;
  const double coverh=ckms/100.;
	
  bin1=(int)arg[0];
  bin2=(int)arg[1];
  bin3=(int)arg[2];

  if (a >= 1.0) error("a>=1 in int_for_b_2x");
  s1 = sglobi;
  s2 = sglobj;
  s3 = sglobk;

  asqr = a*a; 
  kai  = f_K(w(a))*coverh;//kai in h-1Mpc
  f1   = s1/kai;
  f2   = s2/kai;
  f3   = s3/kai;

  hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));
  res = g_sourcex1(a,bin1)*g_sourcex2(a,bin2)*g_sourcex3(a,bin3)/(asqr*asqr*a)/w(a)/hoverh0;//g_sourcex

  if (res>1.e-15) {   //skip if 0 anyway
    k_NL_bs1 = k_NL_bs(a);
    res *= B_3dbs(a,f1,f2,f3,k_NL_bs1)/(ps_ch_rescale*ps_ch_rescale);
  }
  return res;
}


double B_2x_GGI(double s1, double s2, double s3, int bin1, int bin2, int bin3, int iatype)	//tomo B_2_GGI
{
	double arg[4], integral;

	arg[3]=iatype+0.1;

	assert(s1>=s_min);
	assert(s1<=s_max);
	assert(s2>=s_min);
	assert(s2<=s_max);
	assert(s3>=s_min);
	assert(s3<=s_max);


	arg[0]=bin1+0.1;
	arg[1]=bin2+0.1;
	arg[2]=bin3+0.1;

	sglobi = s1;
	sglobj = s2;
	sglobk = s3;
	//printf("begin integrating for zbin: %d-%d-%d; l-bin: %f-%f-%f\n",bin1,bin2,bin3,sglobi,sglobj,sglobk);
	integral=int_GSL_integrate_qag(int_for_b_2x_GGI,arg,a_min,1.0,NULL,2048);
	//printf("end integrating for zbin: %d-%d-%d; l-bin: %f-%f-%f\n",bin1,bin2,bin3,sglobi,sglobj,sglobk);
	//printf("integral=%g\n",integral);
	
	return (9./4.*cosmopar->omm*cosmopar->omm*integral);
}


double int_for_b_2x_GGI(double a, void *args)
{
	double hoverh0, asqr, s1,s2,s3, kai, f1,f2,f3, res;
	int bin1,bin2,bin3;

	double *arg	=	(double*)args;
	const double coverh=ckms/100.;
	
	bin1=(int)arg[0];
	bin2=(int)arg[1];
	bin3=(int)arg[2];
	int iatype=(int)arg[3];

	if (a >= 1.0) error("a>=1 in int_for_b_2x_GGI");
	s1       = sglobi;
	s2       = sglobj;
	s3       = sglobk;

	asqr    = a*a;
	kai     = f_K(w(a))*coverh;//kai in h-1Mpc
	f1       = s1/kai;
	f2       = s2/kai;
	f3       = s3/kai;

	//hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));

	res = g_sourcex1(a,bin1)*g_sourcex2(a,bin2)*probx3(1./a-1.,bin3);

	res = res/asqr/pow(w(a),2);
	res *= B_3d_GGI(a,f1,f2,f3,iatype)/(ps_ch_rescale*ps_ch_rescale);

	return res;
}


double B_3d_GGI(double a, double k1, double k2, double k3,int iatype)
////k in hMpc-1, B_3dbs in (h-1Mpc)^6
//survey/free parameters: A,a_med,k_ref,r,s
{
  double res;

  if (!iatype) {   //power law model
    double k_NL_bs_med,s,r,k_ref,a_med;
	
    s = cosmopar->s_ia;
    r = cosmopar->r_ia;
    k_ref = cosmopar->k_ref;
    a_med = bg_zdistrpar->amed;
    k_NL_bs_med = k_NL_bs(a_med);

    res = cosmopar->A_ia*(-1.)*B_3dbs(a_med,k_ref,k_ref,k_ref,k_NL_bs_med);
    res = res*(pow((k1*k2/k_ref/k_ref),(s-2.))+pow((k2*k3/k_ref/k_ref),(s-2.))+pow((k3*k1/k_ref/k_ref),(s-2.)));
    res = res*pow(a,(-r-2.));
  }
  else {          //linear alignment model
    double ang_k12, ang_k23, ang_k31, q1, q2, q3,k_NL_bs1,iafac;
    double n_eff1, n_eff2, n_eff3;
    double F2_eff123, F2_eff231, F2_eff312;
    double func_a1, func_a2, func_a3; 
    double func_b1, func_b2, func_b3; 
    double func_c1, func_c2, func_c3; 
    double P_NL1, P_NL2, P_NL3;
    const double coverh=ckms/100.;

    n_eff1 = n_eff_bs(k1);
    n_eff2 = n_eff_bs(k2);
    n_eff3 = n_eff_bs(k3);
    
    k_NL_bs1=k_NL_bs(a);
    q1 = k1/k_NL_bs1;	
    q2 = k2/k_NL_bs1;
    q3 = k3/k_NL_bs1;
    
    ang_k12 = -1./2.*(k1/k2+k2/k1-k3*k3/k1/k2);
    ang_k23 = -1./2.*(k2/k3+k3/k2-k1*k1/k2/k3);
    ang_k31 = -1./2.*(k3/k1+k1/k3-k2*k2/k3/k1);

    func_a1 = func_a(a, k1, n_eff1, q1);
    func_a2 = func_a(a, k2, n_eff2, q2);
    func_a3 = func_a(a, k3, n_eff3, q3);

    func_b1 = func_b(a, k1, n_eff1, q1);
    func_b2 = func_b(a, k2, n_eff2, q2);
    func_b3 = func_b(a, k3, n_eff3, q3);

    func_c1 = func_c(a, k1, n_eff1, q1);
    func_c2 = func_c(a, k2, n_eff2, q2);
    func_c3 = func_c(a, k3, n_eff3, q3);

    F2_eff123 = 5./7.*func_a1*func_a2 +1./2.*(k1/k2+k2/k1)*ang_k12*func_b1*func_b2 +2./7.*ang_k12*ang_k12*func_c1*func_c2;

    F2_eff231 = 5./7.*func_a2*func_a3+1./2.*(k2/k3+k3/k2)*ang_k23*func_b2*func_b3 +2./7.*ang_k23*ang_k23*func_c2*func_c3;

    F2_eff312 = 5./7.*func_a3*func_a1+1./2.*(k3/k1+k1/k3)*ang_k31*func_b3*func_b1 +2./7.*ang_k31*ang_k31*func_c3*func_c1;

    P_NL1 = P_NL(a, k1*coverh)*ps_ch_rescale;// P_NL1: dimension length^3, in units of Mpc^3
    P_NL2 = P_NL(a, k2*coverh)*ps_ch_rescale;
    P_NL3 = P_NL(a, k3*coverh)*ps_ch_rescale;

    iafac=C1rhocr*cosmopar->omm/(a*a)/(growfac(a,1.,1.)/growfac(1.,1.,1.));

    res = 2.*iafac*(iafac*F2_eff123*P_NL1*P_NL2 + F2_eff231*P_NL2*P_NL3 + F2_eff312*P_NL3*P_NL1);
  }

  return res;
}



double B_2x_GII(double s1, double s2, double s3, int bin1, int bin2, int bin3, int iatype)	//tomo B_2_GII
{
  double arg[4], integral;

  arg[3]=iatype+0.1;

  assert(s1>=s_min);
  assert(s1<=s_max);
  assert(s2>=s_min);
  assert(s2<=s_max);
  assert(s3>=s_min);
  assert(s3<=s_max);

  
  arg[0]=bin1+0.1;
  arg[1]=bin2+0.1;
  arg[2]=bin3+0.1;
  
  sglobi = s1;
  sglobj = s2;
  sglobk = s3;
  //printf("begin integrating for zbin: %d-%d-%d; l-bin: %f-%f-%f\n",bin1,bin2,bin3,sglobi,sglobj,sglobk);
  integral=int_GSL_integrate_qag(int_for_b_2x_GII,arg,a_min,1.0,NULL,2048);
  //printf("end integrating for zbin: %d-%d-%d; l-bin: %f-%f-%f\n",bin1,bin2,bin3,sglobi,sglobj,sglobk);
  //printf("integral=%g\n",integral);

  return (3./2.*cosmopar->omm*integral);
}


double int_for_b_2x_GII(double a, void *args)
{
	double hoverh0, asqr, s1,s2,s3, kai, f1,f2,f3, res;
	int bin1,bin2,bin3;

	double *arg	=	(double*)args;
	const double coverh=ckms/100.;

	bin1=(int)arg[0];
	bin2=(int)arg[1];
	bin3=(int)arg[2];
	int iatype=(int)arg[3];

	if (a >= 1.0) error("a>=1 in int_for_b_2x_GII");
	s1       = sglobi;
	s2       = sglobj;
	s3       = sglobk;

	asqr    = a*a;
	kai     = f_K(w(a))*coverh;//kai in h-1Mpc
	f1       = s1/kai;
	f2       = s2/kai;
	f3       = s3/kai;

	hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));

	res = g_sourcex1(a,bin1)*probx2(1./a-1.,bin2)*probx3(1./a-1.,bin3);
	res = res*a*hoverh0/pow(w(a),3);

	res *= B_3d_GII(a,f1,f2,f3,iatype)/(ps_ch_rescale*ps_ch_rescale);

	return res;
}


double B_3d_GII(double a, double k1, double k2, double k3, int iatype)
////k in hMpc-1, B_3dbs in (h-1Mpc)^6
//survey;free parameters: a_med;k_ref,r,s,A
{
  double res;

  if (!iatype) {  //power law model
    double k_NL_bs_med,s,r,k_ref,a_med;
	
    s = cosmopar->s_ia;
    r = cosmopar->r_ia;
    k_ref = cosmopar->k_ref;
    a_med = bg_zdistrpar->amed;
    k_NL_bs_med = k_NL_bs(a_med);
    
    res = cosmopar->A_ia*(1./3.)*B_3dbs(a_med,k_ref,k_ref,k_ref,k_NL_bs_med);
    res = res*(pow((k1*k2/k_ref/k_ref),(s-2.))+pow((k2*k3/k_ref/k_ref),(s-2.))+pow((k3*k1/k_ref/k_ref),(s-2.)));
    res = res*pow(a,(-r-2.));
  }
  else {         //linear alignment model
    double ang_k12, ang_k23, ang_k31, q1, q2, q3,k_NL_bs1,iafac;
    double n_eff1, n_eff2, n_eff3;
    double F2_eff123, F2_eff231, F2_eff312;
    double func_a1, func_a2, func_a3; 
    double func_b1, func_b2, func_b3; 
    double func_c1, func_c2, func_c3; 
    double P_NL1, P_NL2, P_NL3;
    const double coverh=ckms/100.;

    n_eff1 = n_eff_bs(k1);
    n_eff2 = n_eff_bs(k2);
    n_eff3 = n_eff_bs(k3);
    
    k_NL_bs1=k_NL_bs(a);
    q1 = k1/k_NL_bs1;	
    q2 = k2/k_NL_bs1;
    q3 = k3/k_NL_bs1;
    
    ang_k12 = -1./2.*(k1/k2+k2/k1-k3*k3/k1/k2);
    ang_k23 = -1./2.*(k2/k3+k3/k2-k1*k1/k2/k3);
    ang_k31 = -1./2.*(k3/k1+k1/k3-k2*k2/k3/k1);

    func_a1 = func_a(a, k1, n_eff1, q1);
    func_a2 = func_a(a, k2, n_eff2, q2);
    func_a3 = func_a(a, k3, n_eff3, q3);

    func_b1 = func_b(a, k1, n_eff1, q1);
    func_b2 = func_b(a, k2, n_eff2, q2);
    func_b3 = func_b(a, k3, n_eff3, q3);

    func_c1 = func_c(a, k1, n_eff1, q1);
    func_c2 = func_c(a, k2, n_eff2, q2);
    func_c3 = func_c(a, k3, n_eff3, q3);

    F2_eff123 = 5./7.*func_a1*func_a2 +1./2.*(k1/k2+k2/k1)*ang_k12*func_b1*func_b2 +2./7.*ang_k12*ang_k12*func_c1*func_c2;

    F2_eff231 = 5./7.*func_a2*func_a3+1./2.*(k2/k3+k3/k2)*ang_k23*func_b2*func_b3 +2./7.*ang_k23*ang_k23*func_c2*func_c3;

    F2_eff312 = 5./7.*func_a3*func_a1+1./2.*(k3/k1+k1/k3)*ang_k31*func_b3*func_b1 +2./7.*ang_k31*ang_k31*func_c3*func_c1;

    P_NL1 = P_NL(a, k1*coverh)*ps_ch_rescale;// P_NL1: dimension length^3, in units of Mpc^3
    P_NL2 = P_NL(a, k2*coverh)*ps_ch_rescale;
    P_NL3 = P_NL(a, k3*coverh)*ps_ch_rescale;

    iafac=C1rhocr*cosmopar->omm/(a*a)/(growfac(a,1.,1.)/growfac(1.,1.,1.));

    res = 2.*iafac*iafac*(iafac*F2_eff123*P_NL1*P_NL2 + F2_eff231*P_NL2*P_NL3 + iafac*F2_eff312*P_NL3*P_NL1);
  }

  return res;
}


double B_2x_III(double s1, double s2, double s3, int bin1, int bin2, int bin3, int iatype)	//tomo B_2_III
{
  double arg[4], integral;

  arg[3]=iatype+0.1;

  assert(s1>=s_min);
  assert(s1<=s_max);
  assert(s2>=s_min);
  assert(s2<=s_max);
  assert(s3>=s_min);
  assert(s3<=s_max);

  
  arg[0]=bin1+0.1;
  arg[1]=bin2+0.1;
  arg[2]=bin3+0.1;
  
  sglobi = s1;
  sglobj = s2;
  sglobk = s3;
  //printf("begin integrating for zbin: %d-%d-%d; l-bin: %f-%f-%f\n",bin1,bin2,bin3,sglobi,sglobj,sglobk);
  integral=int_GSL_integrate_qag(int_for_b_2x_III,arg,a_min,1.0,NULL,2048);
  //printf("end integrating for zbin: %d-%d-%d; l-bin: %f-%f-%f\n",bin1,bin2,bin3,sglobi,sglobj,sglobk);
  //printf("integral=%g\n",integral);

  return (integral);
}


double int_for_b_2x_III(double a, void *args)
{
	double hoverh0, asqr, s1,s2,s3, kai, f1,f2,f3, res;
	int bin1,bin2,bin3;

	double *arg	=	(double*)args;
	const double coverh=ckms/100.;
	
	bin1=(int)arg[0];
	bin2=(int)arg[1];
	bin3=(int)arg[2];
	int iatype=(int)arg[3];

	if (a >= 1.0) error("a>=1 in int_for_b_2x_III");
	s1       = sglobi;
	s2       = sglobj;
	s3       = sglobk;

	asqr    = a*a;
	kai     = f_K(w(a))*coverh;//kai in h-1Mpc
	f1       = s1/kai;
	f2       = s2/kai;
	f3       = s3/kai;

	hoverh0 = sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));

	res = probx1(1./a-1.,bin1)*probx2(1./a-1.,bin2)*probx3(1./a-1.,bin3);
	res *= asqr*asqr*hoverh0*hoverh0/pow(w(a),4);

	res *= B_3d_III(a,f1,f2,f3,iatype)/(ps_ch_rescale*ps_ch_rescale);

	return res;
}


double B_3d_III(double a, double k1, double k2, double k3, int iatype)
////k in hMpc-1, B_3dbs in (h-1Mpc)^6
//survey;free parameters: a_med;k_ref,r,s,A
{
  double res;

  if (!iatype) {  //power law model
    double k_NL_bs_med,s,r,k_ref,a_med;
	
    s = cosmopar->s_ia;
    r = cosmopar->r_ia;
    k_ref = cosmopar->k_ref;
    a_med = bg_zdistrpar->amed;
    k_NL_bs_med = k_NL_bs(a_med);
    
    res = cosmopar->A_ia*(1./3.)*B_3dbs(a_med,k_ref,k_ref,k_ref,k_NL_bs_med);
    res = res*(pow((k1*k2/k_ref/k_ref),(s-2.))+pow((k2*k3/k_ref/k_ref),(s-2.))+pow((k3*k1/k_ref/k_ref),(s-2.)));
    res = res*pow(a,(-r-2.));
  }
  else {         //linear alignment model
    double iafac;

    iafac=C1rhocr*cosmopar->omm/(a*a)/(growfac(a,1.,1.)/growfac(1.,1.,1.));
    res=iafac*iafac*iafac*iafac*B_3dbs(a,k1,k2,k3,k_NL_bs(a));
  }

  return res;
}




// Calculation of Map^3 and Map^2, added by Benjamin
#define nlintbin 40     // global definitions for map3-integration
#define maxint 20.      // integrate fouriermapfilter from minint to maxint
#define minint 1.e-5


double fouriermapfilter(double x,int filter)
{
  double x2=x*x;
  if (filter==0) return(12./PI*gsl_sf_bessel_Jn(4,x)/x2);
  else if (filter==1) return(x2*exp(-x2/2.)/(4.*PI));
  else {
    printf("Error in routine fouriermapfilter: unknown value for filter type: %i\n",filter);
    exit(-1);
  }
}

double map3_Ia(double theta,double l1,double l2)
{
  double tsq=theta*theta;
  double lsq=l1*l1+l2*l2;
  double arg=l1*l2*tsq;

  return(0.25*exp(-0.5*lsq*tsq)*tsq*(lsq*gsl_sf_bessel_In(0,arg)-2.*l1*l2*gsl_sf_bessel_In(1,arg)));
}

double map3_Ib(double theta,double l1,double l2)
{
  double tsq=theta*theta;
  double lsq=l1*l1+l2*l2;
  double arg=l1*l2*tsq;

  return(exp(-0.5*lsq*tsq)*((0.5-0.25*lsq*tsq)*gsl_sf_bessel_In(1,arg)+0.5*arg*gsl_sf_bessel_In(2,arg)));
}

double map3_Ic(double theta,double l1,double l2)
{
  double tsq=theta*theta;
  double lsq=l1*l1+l2*l2;
  double arg=l1*l2*tsq;

  return(exp(-0.5*lsq*tsq)/(4.*PI*l1*l2)*(PI*lsq*gsl_sf_bessel_In(1,arg)+l1*l2*((-6.*PI+PI*lsq*tsq)*gsl_sf_bessel_In(2,arg)-2.*PI*arg*gsl_sf_bessel_In(3,arg))));
}


double int_map3_S(double a, void *args)
{
  double *par=(double*)args;  //arg[0]=l1,arg[1]=l2,arg[2]=signaltype,arg[3]=fitformulaterm
  double hoverh0,asqr,s1,s2,kai,k1,k2,res=0.,k_NL_bs1,g,p,chi,iafac,ffterm=1.;
  double q1,q2;
  double n_eff1, n_eff2;
  double P1, P2;
  const double coverh=ckms/100.;

  if (a >= 1.0) error("a>=1 in int_for_map3_Sx");
  asqr=a*a;

  s1=par[0];
  s2=par[1];
  int iatype=(int)par[2];
  int fftype=(int)par[3];

  kai=f_K(w(a))*coverh;//kai in h-1Mpc
  k1=s1/kai;
  k2=s2/kai;

  hoverh0=sqrt(cosmopar->omm/(a*asqr) + (1.-cosmopar->omm-cosmopar->omv)/asqr + OMV_MACRO(a));
  g=g_source(a);
  p=prob(1./a-1.);
  chi=w(a);
  iafac=C1rhocr*cosmopar->omm/asqr/(growfac(a,1.,1.)/growfac(1.,1.,1.));

  if (iatype==0) res=g*g*g/(asqr*asqr*a)/chi/hoverh0;
  else if (iatype==1) res=p*g*g/asqr/(chi*chi)*iafac*(2.+iafac);
  else if (iatype==2) res=p*p*g*a*hoverh0/(chi*chi*chi)*iafac*iafac*(1.+2.*iafac);
  else res=p*p*p*asqr*asqr*hoverh0*hoverh0/(chi*chi*chi*chi)*iafac*iafac*iafac*iafac;

  k_NL_bs1 = k_NL_bs(a);

  n_eff1 = n_eff_bs(k1);
  n_eff2 = n_eff_bs(k2);

  q1 = k1/k_NL_bs1;	
  q2 = k2/k_NL_bs1;

  if (cosmopar->iamethod==2) {
    printf("Error in routine int_map3_S: IA halo model not implemented in this function.\n");
    exit(-1);
  }

  if ((fftype==0)&&((iatype==0)||(cosmopar->iamethod==1))) ffterm=func_a(a, k1, n_eff1, q1)*func_a(a, k2, n_eff2, q2);
  else if ((fftype==1)&&((iatype==0)||(cosmopar->iamethod==1))) ffterm=func_b(a, k1, n_eff1, q1)*func_b(a, k2, n_eff2, q2)*(k1/k2+k2/k1);
  else if ((fftype==2)&&((iatype==0)||(cosmopar->iamethod==1))) ffterm=func_c(a, k1, n_eff1, q1)*func_c(a, k2, n_eff2, q2);

  if ((iatype==0)||(cosmopar->iamethod==1)) {
    P1 = P_NL(a, k1*coverh);
    P2 = P_NL(a, k2*coverh);  // Power spectrum: *pow(3000., 3) to get units of Mpc^3
  }
  else {
    P1 = P_L(a, k1*coverh);
    P2 = P_L(a, k2*coverh);
  }
  return(ffterm*P1*P2*res);
}


double int_map3_l2(double l2, void *args) 
{
  double *arg=(double*)args;  //arg[0]=theta,arg[1]=l1,arg[2]=signaltype,arg[3]=filtertype
  double par[4];
  par[2]=arg[2];
  double lintmax=log(maxint/arg[0]);
  double lintmin=log(minint/arg[0]);
  double dl=(lintmax-lintmin)/(nlintbin-1.);

  static double THETA=-42.;
  static double ***table;

  if (arg[0]!=THETA) {

    int i,j;
    table=malloc(3*sizeof(double *));
    for(i=0;i<3;i++) {
      table[i]=malloc(nlintbin*sizeof(double *));
      for(j=0;j<nlintbin;j++) {
	table[i][j]=malloc(nlintbin*sizeof(double));
      }
    }

    for(i=0;i<nlintbin;i++) {
      par[0]=exp(lintmin+i*dl);
      for(j=0;j<nlintbin;j++) {
	par[1]=exp(lintmin+j*dl);

	par[3]=0.1;
	table[0][i][j]=int_GSL_integrate_qag(int_map3_S,par,a_min,1.0,NULL,2048);  //a

	par[3]=1.1;
	table[1][i][j]=int_GSL_integrate_qag(int_map3_S,par,a_min,1.0,NULL,2048);  //b

	par[3]=2.1;
	table[2][i][j]=int_GSL_integrate_qag(int_map3_S,par,a_min,1.0,NULL,2048);  //c

      }
    }
    THETA=arg[0];
  }
  
  return(fouriermapfilter(arg[0]*l2,(int)arg[3])*l2*(5./7.*map3_Ia(arg[0],arg[1],l2)*interpol2d(table[0],nlintbin,lintmin,lintmax,dl,log(arg[1]),nlintbin,lintmin,lintmax,dl,log(l2),0.,0.) + 0.5*map3_Ib(arg[0],arg[1],l2)*interpol2d(table[1],nlintbin,lintmin,lintmax,dl,log(arg[1]),nlintbin,lintmin,lintmax,dl,log(l2),0.,0.) + 2./7.*map3_Ic(arg[0],arg[1],l2)*interpol2d(table[2],nlintbin,lintmin,lintmax,dl,log(arg[1]),nlintbin,lintmin,lintmax,dl,log(l2),0.,0.)));
}

double int_map3_l1(double l1, void *args) 
{
  double *arg=(double*)args;  //arg[0]=theta,arg[1]=signaltype,arg[2]=filtertype
  double par[4];
  par[0]=arg[0];
  par[2]=arg[1];
  par[3]=arg[2];
  double lintmax=log(maxint/arg[0]);
  double lintmin=log(minint/arg[0]);
  double dl=(lintmax-lintmin)/(nlintbin-1.);

  static double THETA=-42.;
  static double table[nlintbin];

  if (arg[0]!=THETA) {

    int i;
    for(i=0;i<nlintbin;i++) {
      par[1]=exp(lintmin+i*dl);
      table[i]=int_GSL_integrate_qag(int_map3_l2,par,exp(lintmin),exp(lintmax),NULL,2048);
    }
    THETA=arg[0];
  }

  return(fouriermapfilter(arg[0]*l1,(int)arg[2])*l1*interpol(table,nlintbin,lintmin,lintmax,dl,log(l1),0.,0.));
}


double map3(double theta,int signaltype,int filtertype) 
{
  double par[3]={theta,signaltype+.1,filtertype+.1};
  double prefac=12.*pow(1.5*cosmopar->omm,3.-1.*signaltype);  //*3 for 3 permutations of bispectrum; *2 due to prefactor in front of F(k1,k2); *2 due to integration over phi up to Pi instead of 2*Pi 

  return(prefac*int_GSL_integrate_qag(int_map3_l1,par,minint/theta,maxint/theta,NULL,2048));
}


double int_map2(double l, void *args) 
{
  double *arg=(double*)args;  //arg[0]=theta,arg[1]=bin1,arg[2]=bin2,arg[3]=signaltype,arg[4]=filtertype
  int bin1=(int)arg[1];
  int bin2=(int)arg[2];
  int type=(int)arg[3];
  int filtertype=(int)arg[4];
  double ps=0.;

  double filter=fouriermapfilter(l*arg[0],filtertype);
  if (type==0) ps=P_2x(l,bin1,bin2);
  else if (type==1) ps=P_2iax(l,bin1,bin2);
  else ps=P_2IIx(l,bin1,bin2);

  return(l*ps*filter*filter);
}


 double map2(double theta,int bin1,int bin2,int signaltype,int filtertype) 
{
  double par[5]={theta,bin1+.1,bin2+.1,signaltype+.1,filtertype+.1};

  if (bg_zdistrpar->photo_sig==0.0 && bin1!=bin2) return 0.0;
  else return(2*PI*int_GSL_integrate_qag(int_map2,par,minint/theta,maxint/theta,NULL,2048));
}



// procedures for computing GI correlation functions
double int_for_gixi(double k, void *args)
{
  double *arg = (double*)args;  //0: x_perp; 1: a
  return(k*gsl_sf_bessel_Jn(2,k*arg[0])*P_GI(arg[1],k));
}


double gixi(double x, double z)
{

  #define STOP 1e-4         // relative contribution per integration interval
  #define INTERV 100        // Scaling of integration range

  static double OMEGA_M = -42.;
  static double OMEGA_V = -42.;
  static double W0      = -42.;
  static double WA      = -42.;
  static double N_SPEC  = -42.;
  static double GAMMA   = -42.;
  static double BETA_P  = -42.;
  static double ZMAX    = -42.;
  static double ZMIN    = -42.;
  static double Z0      = -42.;
  static double FCAT    = -42.;
  static double SIG     = -42.;
  static double DELTAZ  = -42.;
  static double NONLINEAR = -42.;
  static double SIGMA_8 = -42.;
  static double Z = -42.;
  static int PSNLTYPE=-42;
  static int EH=-42;
  static double table[N_s];
  static double dx = 0.0;
  double logxmin, logxmax;
  int    i,j;
  double arg[2], integral,diff,res;

  const double xmin=1.e-1;
  const double xmax=1.e2;
  logxmin = log(xmin);
  logxmax = log(xmax);

  if (OMEGA_M != cosmopar->omm || OMEGA_V != cosmopar->omv || W0 != cosmopar->w0 || WA != cosmopar->wa || N_SPEC != cosmopar->n || GAMMA != cosmopar->Gamma || BETA_P != bg_zdistrpar->beta || Z0 != bg_zdistrpar->z0 || FCAT != bg_zdistrpar->photo_fcat || SIG != bg_zdistrpar->photo_sig || DELTAZ != bg_zdistrpar->photo_deltaz || NONLINEAR != cosmopar->nonlinear || SIGMA_8 != cosmopar->sigma8 || ZMAX != bg_zdistrpar->zmax || ZMIN != bg_zdistrpar->zmin || PSNLTYPE != cosmopar->psnltype || EH != cosmopar->transfer_EH || Z != z)	
  {
    if(N_s==1) dx=0.; 
    else dx = (logxmax - logxmin)/(1.*N_s - 1.);

    arg[1]=1./(1.+z);

    for (i=0;i<N_s;i++) {
      arg[0] = exp(logxmin+i*dx);
      diff=INTERV/(arg[0]+1.);
      integral=0.0;
      j=0;
      do { 
	res=int_GSL_integrate_qag(int_for_gixi,arg,j*diff,(j+1)*diff,NULL,2048);	
	integral+=res;
	printf("%i  %g\n",j,fabs(res/integral));
	j++;
      } while (fabs(res/integral)>STOP);
      table[i] = (-1.)*integral/(2.*PI);
    }
    OMEGA_M =   cosmopar->omm;
    OMEGA_V =   cosmopar->omv;
    W0      =   cosmopar->w0;
    WA      =   cosmopar->wa;
    N_SPEC  =   cosmopar->n;
    GAMMA   =   cosmopar->Gamma;
    BETA_P  =   bg_zdistrpar->beta;
    Z0      =   bg_zdistrpar->z0;
    FCAT    =   bg_zdistrpar->photo_fcat;
    SIG     =   bg_zdistrpar->photo_sig;
    DELTAZ  =   bg_zdistrpar->photo_deltaz;
    NONLINEAR = cosmopar->nonlinear;
    SIGMA_8 =   cosmopar->sigma8;
    ZMAX    =   bg_zdistrpar->zmax;
    ZMIN    =   bg_zdistrpar->zmin;
    PSNLTYPE   =	cosmopar->psnltype;
    EH	    =	cosmopar->transfer_EH;
    Z       =   z;
  }

  return(interpol(table, N_s, logxmin, logxmax, dx, log(x), 0., 0.));
}

