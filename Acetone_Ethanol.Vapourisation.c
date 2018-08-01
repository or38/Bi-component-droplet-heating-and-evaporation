/**********************************************************************
UDF for defining the heat and mass transport for Bi-component Droplet vaporization
Has been tested with ANSYS Fluent 17.2

Sazhin, S.S., et al., A simplified model for bi-component droplet heating and evaporation. International Journal of Heat and Mass Transfer, 2010. 53(21-22): p. 4495-4505.


If you use our software, please cite the following:
1. O. Rybdylova, M. Al Qubeissi, M. Braun, C. Crua, J. Manin, L.M. Pickett, G. de Sercey, E.M. Sazhina, S.S. Sazhin, M. Heikal, A model for droplet heating and its implementation into ANSYS Fluent, International Communications in Heat and Mass Transfer, Volume 76, 2016, Pages 265-270, ISSN 0735-1933,
https://doi.org/10.1016/j.icheatmasstransfer.2016.05.032.

2. O. Rybdylova, L. Poulton, M. Al Qubeissi, A.E. Elwardany, C. Crua, T. Khan, S.S. Sazhin, A model for multi-component droplet heating and evaporation and its implementation into ANSYS Fluent, International Communications in Heat and Mass Transfer, Volume 90, 2018, Pages 29-33, ISSN 0735-1933,
https://doi.org/10.1016/j.icheatmasstransfer.2017.10.018.

Copyright (C) 2018 Oyuna Rybdylova, Luke Poulton - All Rights Reserved
You may use, distribute and modify this code under the terms of the MIT license

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
***********************************************************************/
#include "udf.h"
#include <time.h>

#define BM_MAX 1.E20
#define BM_MIN -0.99999
#define ACCURACY 1.e-6
#define PI 3.1415926535897932384626433832795
#define R_AIR 287.01625988193461525183829875375
#define R_UNIVERSAL 8.3144622e3
#define N_Lambda 200 // number of terms in the series
#define N_Lambda_y 10 // number of terms in the series for mass fraction
#define N_INT 500 // number of layers inside a droplet for temperature and each component
#define Delta_R 0.002 // = 1/N_INT
#define N_spieces 2
#define Ind(i,j) i+N_spieces*j
#define N_T_average 1503//303
#define N_x_0 1506//306
#define N_components 2
#define N_LJ_l 1512//312
#define N_mdot 1516//316

/* convection diffusion controlled vaporisation model as implemented into Fluent
   p    ... tracked particle struct
   Cp   ... particle heat capacity
   hgas ... enthalpy of formation for gas species
   hvap ... vaporization enthalpy
   cvap_surf ... molar concentration of vapor at surface
   Z    ... compressibility factor
   dydt ... temperature / component mass source terms for particle 
            dydt[0] = temperature equation, 
            dydt[1+ns] = mass equation of particle species ns
   dzdt ... temperature / species mass source terms for gas phase
            dzdt->energy = convective heat transfer in Watt
            dzdt->species[n_gas] = mass source for gas species n_gas in kg/s
 */

//N_INT + 1 - Temperature, (N_INT + 1 )*(N_components) - mass fractions + 1 - T_average + N_component - Y_i_average + N_component x_svi, N_component - gamma_i , N_component - e_i vap rate, N_component Lennard-Jones lengths, N_component - LJ energies, m_dot, N_componnents M_idot , N_components x_sli, N_components y_vsi//322
/*0-100 Temperature*/
/*101-201 Y_0*/
/*202-302 Y_1*/
/*303 T_average*/
/*304, 305 Y_i_average*/
/*306, 307 X_i_vapour*/
/*308, 309 gamma_i*/
/*310, 311 evaporation rate_i (epsilon_i)*/
/*312, 313 Lennard-Jones length_i*/
/*314, 315 Lennard-Jones energy*/
/*316 total evaporation rate*/
/*317, 318 m_dot_i*/
/*319, 320 x_ls*/
/*321, 322 Y_vap_i*/
/*323 in total*/

/*0-500 Temperature*/
/*501-1001 Y_0*/
/*1002-1502 Y_1*/
/*1503 T_average*/
/*1504, 1505 Y_i_average*/
/*1506, 307 X_i_vapour*/
/*1508, 309 gamma_i*/
/*1510, 311 evaporation rate_i (epsilon_i)*/
/*1512, 313 Lennard-Jones length_i*/
/*1514, 315 Lennard-Jones energy*/
/*1516 total evaporation rate*/
/*1517, 318 m_dot_i*/
/*1519, 320 x_ls*/
/*1521, 322 Y_vap_i*/
/*1523 in total*/

int Lambda(real h_0, real lambda[])
{
	//FILE * fout;
	int i, j, indexi;
	int Max_iterations = 200;
	real lambda_left, lambda_right, f_left, f_right, lambda_mid, f_mid;
	real conv_crit = 1.e-13;
	real step = 1.e-7;

	for (i = 0; i < N_Lambda; i++) lambda[i] = -1.0;

	for (i = 0; i < N_Lambda; i++)
	{
		indexi = 0;
        j = i;
        if(fabs(h_0)>1.0) j = i + 1;
        //Message("j = %d\n", j);
        //Message("h_0 = %f\tabs(h_0) = %f\n", h_0,fabs(h_0));

        lambda_left = ((real)(j))*PI + step;
		lambda_right = (((real)(j + 1)) - 0.5)*PI - step;

		if (h_0 > 0.0)
		{
			lambda_left += 0.5*PI;
			lambda_right += PI*0.5;
		}

		f_left = lambda_left*cos(lambda_left) + h_0*sin(lambda_left);
		f_right = lambda_right*cos(lambda_right) + h_0*sin(lambda_right);
        //Message("Lambda left = %f\nLambda right = %f\n f left = %e\nf right = %e\n------------\n", lambda_left, lambda_right, f_left, f_right);
		if (f_left*f_right < 0.0)
		{
			while ((lambda_right - lambda_left > conv_crit) && (indexi<Max_iterations))
			{
				lambda_mid = (lambda_left + lambda_right)*0.5;
				f_mid = lambda_mid*cos(lambda_mid) + h_0*sin(lambda_mid);
				if (f_left*f_mid < 0.0)
				{
					lambda_right = lambda_mid;
					f_right = lambda_right*cos(lambda_right) + h_0*sin(lambda_right);
				}
				else
				{
					lambda_left = lambda_mid;
					f_left = lambda_left*cos(lambda_left) + h_0*sin(lambda_left);
				}
				//Message("Lambda left = %f\nLambda right = %f\n f left = %e\nf right = %e\n------------\n", lambda_left, lambda_right, f_left, f_right);
				indexi++;
			}
			lambda[i] = lambda_left;
			//Message("Lambda[%d]: Number of iterations is %d\n", i,indexi);
		}
   // Message("Lambda[%d] = %f\t", i,lambda[i]);    
	}
    //Message("\n");
	//	Message("Lambdas are calculated\n");

	//	fout = fopen("lambda.txt", "w");
	//	for (i = 0; i < N_Lambda; i++)
	//		fprintf(fout, "%d\t%20.19f\t%e\n", i + 1, lambda[i], lambda[i] * cos(lambda[i]) + h_0*sin(lambda[i]));
	//	fclose(fout);
	//	Message("Lambdas are printed in lambda.txt\n");
	return 0;
}
int Lambda_y(real h_0, real lambda[])
{
	//FILE * fout;
	int i, j, indexi;
	int Max_iterations = 200;
	real lambda_left, lambda_right, f_left, f_right, lambda_mid, f_mid;
	real conv_crit = 1.e-13;
	real step = 1.e-7;

	for (i = 0; i < N_Lambda_y; i++) lambda[i] = -1.0;

	for (i = 0; i < N_Lambda_y; i++)
	{
		indexi = 0;
		j = i;
		if (fabs(h_0)>1.0) j = i + 1;
		//Message("j = %d\n", j);
		//Message("h_0 = %f\tabs(h_0) = %f\n", h_0,fabs(h_0));

		lambda_left = ((real)(j))*PI + step;
		lambda_right = (((real)(j + 1)) - 0.5)*PI - step;

		if (h_0 > 0.0)
		{
			lambda_left += 0.5*PI;
			lambda_right += PI*0.5;
		}

		f_left = lambda_left*cos(lambda_left) + h_0*sin(lambda_left);
		f_right = lambda_right*cos(lambda_right) + h_0*sin(lambda_right);
		//Message("Lambda left = %f\nLambda right = %f\n f left = %e\nf right = %e\n------------\n", lambda_left, lambda_right, f_left, f_right);
		if (f_left*f_right < 0.0)
		{
			while ((lambda_right - lambda_left > conv_crit) && (indexi<Max_iterations))
			{
				lambda_mid = (lambda_left + lambda_right)*0.5;
				f_mid = lambda_mid*cos(lambda_mid) + h_0*sin(lambda_mid);
				if (f_left*f_mid < 0.0)
				{
					lambda_right = lambda_mid;
					f_right = lambda_right*cos(lambda_right) + h_0*sin(lambda_right);
				}
				else
				{
					lambda_left = lambda_mid;
					f_left = lambda_left*cos(lambda_left) + h_0*sin(lambda_left);
				}
				//Message("Lambda left = %f\nLambda right = %f\n f left = %e\nf right = %e\n------------\n", lambda_left, lambda_right, f_left, f_right);
				indexi++;
			}
			lambda[i] = lambda_left;
			//Message("Lambda[%d]: Number of iterations is %d\n", i,indexi);
		}
		Message("Lambda[%d] = %f\t", i,lambda[i]);    
	}
	//Message("\n");
	//	Message("Lambdas are calculated\n");

	//	fout = fopen("lambda.txt", "w");
	//	for (i = 0; i < N_Lambda; i++)
	//		fprintf(fout, "%d\t%20.19f\t%e\n", i + 1, lambda[i], lambda[i] * cos(lambda[i]) + h_0*sin(lambda[i]));
	//	fclose(fout);
	//	Message("Lambdas are printed in lambda.txt\n");
	return 0;
}
real Lambda_m0(double h_0)
{
	//FILE * fout;
	int i;
	real lambda_left, lambda_right, f_left, f_right, lambda_mid, f_mid;
	real conv_crit = 1.e-13;
	real step = 1.e-7;
    real lambda_m0 = -1.0;

	lambda_left = step;
	lambda_right = 100.0;
	
	f_left = lambda_left*cosh(lambda_left) + h_0*sinh(lambda_left);
	f_right = lambda_right*cosh(lambda_right) + h_0*sinh(lambda_right);
    //Message("Lambda left = %f\nLambda right = %f\n f left = %e\nf right = %e\n------------\n", lambda_left, lambda_right, f_left, f_right);
	if (f_left*f_right < 0.0)
	{
		while (lambda_right - lambda_left > conv_crit)
		{
			lambda_mid = (lambda_left + lambda_right)*0.5;
			f_mid = lambda_mid*cosh(lambda_mid) + h_0*sinh(lambda_mid);
			if (f_left*f_mid < 0.0)
			{
				lambda_right = lambda_mid;
				f_right = lambda_right*cosh(lambda_right) + h_0*sinh(lambda_right);
			}
			else
			{
				lambda_left = lambda_mid;
				f_left = lambda_left*cosh(lambda_left) + h_0*sinh(lambda_left);
			}
			//Message("Lambda left = %f\nLambda right = %f\n f left = %e\nf right = %e\n------------\n", lambda_left, lambda_right, f_left, f_right);
		}
		lambda_m0 = lambda_left;
	}
   // Message("l_m0 = %f\n",lambda_m0);
	return lambda_m0;
}

real Vapour_Density(Tracked_Particle *p, real T, real P)
{
	real y1, y2, y3, r1, r2, r3;
    cphase_state_t *c = &(p->cphase);
	
	y1 = (2.0*P_USER_REAL(p, N_mdot + 5)+c->yi[1])/3.0;
	y2 = (2.0*P_USER_REAL(p, N_mdot + 6)+c->yi[0])/3.0;
	y3 = 1.0 - y1 - y2; //air
	r1 = P * solver_par.molWeight[1] / R_UNIVERSAL / T;
	r2 = P * solver_par.molWeight[0] / R_UNIVERSAL / T;
	r3 = P / R_AIR / T;

	//return 1.0/((y1/r1)+(y2/r2)+(y3/r3));
	return 0.0000134*T*T - 0.0120250*T + 3.5647000;
}
int gamma_activity(Tracked_Particle *p, real T)
{
	double a, b, c;
	real x1, x2, g1,g2;
	Material *cond_mix = P_MATERIAL(p);

	a = 546.3 / T - 0.9897;
	b = 543.3 / T - 0.9483;
	c = 15.63 / T + 0.0759;

	if (cond_mix->component[0]->name[0] == 'e')
	{
		x1 = P_USER_REAL(p, N_mdot + 3);
		x2 = P_USER_REAL(p, N_mdot + 3 + 1);
	}
	else{
		x1 = P_USER_REAL(p, N_mdot + 3 + 1);
		x2 = P_USER_REAL(p, N_mdot + 3);
	}
	g1 = exp(x2*x2*(a + 2.0*(b - a - c)*x1 + 3.0*c*x1*x1)); //eth
	g2 = exp(x1*x1*(b + 2.0*(a - b - c)*x2 + 3.0*c*x2*x2)); //ace
	//Message("g1=%e g2=%e\n", g1, g2);
	
	P_USER_REAL(p, N_x_0 + N_components) = g1; //Ace
	P_USER_REAL(p, N_x_0 + N_components + 1) = g2; //eth

	if (cond_mix->component[0]->name[0] == 'a'){
		P_USER_REAL(p, N_x_0 + N_components) = g2;
		P_USER_REAL(p, N_x_0 + N_components + 1) = g1;
	}

	return 0;
}
real vapour_diffusivity(Tracked_Particle *p, real T, real P)
{
	real x_tot = 0.e-15;
	real Mv = 0.e-15;
	real Mva = 0.e-15;
	real T_n;
	real Omega_D, molwt[2];
	int ns;
	int nc = TP_N_COMPONENTS(p);
	real sigmava, sigma_v, e_v, e_va, D;

	sigma_v = 0.e-15;
	e_v = 0.e-15;
    molwt[1] = 46.069; 
    molwt[0] = 58.08;
	

	for (ns = 0; ns < nc; ns++)
	{
		Mv += P_USER_REAL(p, N_x_0 + ns)*molwt[ns];
		sigma_v += P_USER_REAL(p, N_x_0 + ns)*P_USER_REAL(p, N_LJ_l + ns);
		e_v += P_USER_REAL(p, N_x_0 + ns)*P_USER_REAL(p, N_LJ_l + N_components + ns);
		x_tot += P_USER_REAL(p, N_x_0 + ns);
	}
	Mv = Mv / x_tot;
	sigma_v = sigma_v / x_tot;
	e_v = e_v / x_tot;

	Mva = 2.0 / (1 / Mv + 1 / 28.967); //air
	Mva = sqrt(Mva);
	sigmava = 0.5*(sigma_v + 3.711); //air
	T_n = T / sqrt(78.6*e_v); //air
	//Omega_D = 1.06036*pow(T_n, -0.1561) + 0.193*exp(0.47635*T_n) + 1.03587*exp(1.52996*T_n) + 1.76474*exp(3.89411*T_n);
	Omega_D = 1.06036 / pow(T_n, 0.1561) + 0.193 / exp(0.47635*T_n) + 1.03587 / exp(1.52996*T_n) + 1.76474 / exp(3.89411*T_n);
    //Message("Mv = %f\tP = %f\tMva = %f\tsigmava = %f\tOmega_D = %f\tT = %f\n", Mv, P, Mva, sigmava,Omega_D, T);
    //Message("D = %e\n", (3.03 - 0.98 / Mva) / P / Mva / sigmava / sigmava / Omega_D*1.e-2*pow(T, 1.5));
    //Message("D1 = %e\n", (3.03 - 0.98 / Mva) / P);
    //Message("D = %e\n", 1/ Mva / sigmava / sigmava / Omega_D*1.e-2*pow(T, 1.5));
	//return (3.03 - 0.98 / Mva) / P / Mva / sigmava / sigmava / Omega_D*1.e-2*pow(T, 1.5);
	return D = (3.03 - 0.98 / Mva) / (P * Mva * (sigmava * sigmava) * Omega_D)*1.e-2*pow(T, 1.5);  //wilke and lee//!!!!!!!!!!!!!!!!!
}

real vapour_heatcap(Tracked_Particle *p, real T)
{
	real y1, y2, y3, c1, c2, c3;
	y1 = P_USER_REAL(p, N_mdot + 5); //Acetone??
	y2 = P_USER_REAL(p, N_mdot + 6);
	y3 = 1.0 - y1 - y2; //air
	c1 = 7.0473444e-6*T*T*T - 9.9229425e-3*T*T + 8.211229*T - 458.00814;  //acetone
	c2 = 469.67 + T*(4.2301 - (1.5571e-3)*T); //ethanol
	//c3 = 1024.9 + T * (T * (5.7628e-4 - 2.4225e-7 * T) - 0.21674);//air
	c3 = (0.0000002*T*T - 0.0000900*T + 1.0160000)*1000.0;
	return y1*c1 + y2*c2 + y3*c3 ;
	//Message("Cp.A=%f Cp.E=%f Cpair=%\n", c1, c2, c3);
}
real Wilke_vap_thermcond(Tracked_Particle *p, real T)
{
	real sum = 0.;
	real phi[9];
	real visc[3], molwt[3], tc[3], x[3];
	real Phi, M_tot;
	int i, j;
	cphase_state_t *c = &(p->cphase);  /* continuous phase struct, caching variables of the cell */
	Thread *t = P_CELL_THREAD(p);
	Material *spi;
	Material *spj;
	Material *gas_mix = THREAD_MATERIAL(t);
	Property *prop;
	real result = 0.e-15;

	M_tot = 0.e-15;
	tc[1] = 0.01143468*pow(T / 273.15, 2.0 - 0.342843365e-6*T*T); //acetone
	tc[0] = 1.8037e-3 + 7.419e-6*T + 1.1536e-7*T*T; //ethanol
	tc[2] = 0.0242; //air

	//visc[1] = -0.1579884447e-5 + 0.3090958433e-7*T - 0.449329e-11*T*T;//acetoneoriginal
	visc[1] = -0.1579884447e-5 + 0.3090958433e-7*T - 0.449329e-11*T*T;
	visc[0] = -0.19757e-6 + 0.29211e-7*T;//ethanol
	visc[2] = 1.7894e-05;//air

	molwt[1] = 58.08;
	molwt[0] = 46.07;
	molwt[2] = 28.967;

	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		M_tot += c->yi[i] / molwt[i];
	}

	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			if (i == j) phi[Ind(i, j)] = 1.0;
			else phi[Ind(i, j)] = (1.0 + sqrt(visc[i] / visc[j] * sqrt(molwt[j] / molwt[i]))) / sqrt(8.0) / (sqrt(1.0 + molwt[i] / molwt[j]));
		}
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		x[i] = c->yi[i] / molwt[i] / M_tot;
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		Phi = 0.e-15;
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			Phi += x[j] * phi[Ind(i, j)];
		}
		result += x[i] * tc[i] / Phi;
	}


	return result;
}
real Wilke_vap_visc(Tracked_Particle *p, real T)
{
	real sum = 0.;
	real phi[9];
	real visc[3], molwt[3], x[3];
	real Phi, M_tot;
	int i, j;
	cphase_state_t *c = &(p->cphase);
	Thread *t = P_CELL_THREAD(p);
	Material *spi;
	Material *spj;
	Material *gas_mix = THREAD_MATERIAL(t);
	Property *prop;
	real result = 0.e-15;

	M_tot = 0.e-15;
	visc[0] = -0.19757e-6 + 0.29211e-7*T; //ethanol
	visc[1] = -0.1579884447e-5 + 0.3090958e-7*T - 0.449329e-11*T*T; //acetone
	visc[2] = 1.85497e-5;

	molwt[0] = 46.07; //ethanol
	molwt[1] = 58.08; //acetone
	molwt[2] = 28.967;
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		//prop = (MATERIAL_PROPERTY(spi));
		//visc[i] = generic_property(c, t, prop, PROP_mu, C_T(c, t));
		//molwt[i] = MATERIAL_PROP(MIXTURE_COMPONENT(gas_mix, i), PROP_mwi);
		M_tot += c->yi[i] / molwt[i];
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			if (i == j) phi[Ind(i, j)] = 1.0;
			else phi[Ind(i, j)] = (1.0 + sqrt(visc[i] / visc[j] * sqrt(molwt[j] / molwt[i]))) / sqrt(8.0) / (sqrt(1.0 + molwt[i] / molwt[j]));
		}
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		x[i] = c->yi[i] / molwt[i] / M_tot;
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		Phi = 0.e-15;
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			Phi += x[j] * phi[Ind(i, j)];
		}
		result += x[i] * visc[i] / Phi;
	}

	return result;

}
real LatentHeat(Tracked_Particle *p, real T)
{
	Material *cond_mix = P_MATERIAL(p);
	double k1, k2;

	if (cond_mix->component[0]->name[0] == 'e')
	{
		k1 = 120.91e3*pow(516.2 - T, 0.38); //eth
		k2 = 489.0 * pow((508.1 - T) / (508.1 - 329.22), 0.38)*1000.0; //ace
	}
	else{
		k2 = 120.91e3*pow(516.2 - T, 0.38);//eth
		k1 = 489.0 * pow((508.1 - T) / (508.1 - 329.22), 0.38)*1000.0; //ace
	}
	return P_USER_REAL(p, N_x_0 + 2 * N_components)*k1 + P_USER_REAL(p, N_x_0 + 2 * N_components + 1)*k2;
}
real liquid_thermal_conductivity(Tracked_Particle *p, real T)
{
	Material *cond_mix = P_MATERIAL(p);
	double k1, k2;

	if (cond_mix->component[0]->name[0] == 'e')
	{
		k1 = 0.61572 - 0.24127e-2*T + 0.31333e-5*T*T;//eth
		k2 = 0.3133614225 - 0.8163e-3*T + 0.1e-5*T*T;//ace
	    
	}
	else{
		k2 = 0.61572 - 0.24127e-2*T + 0.31333e-5*T*T;//eth
		k1 = 0.3133614225 - 0.8163e-3*T + 0.1e-5*T*T;//ace
		
	}
	return P_USER_REAL(p, N_T_average + 1)*k2 + P_USER_REAL(p, N_T_average + 2)*k1 - 0.72*P_USER_REAL(p, N_T_average + 1)*P_USER_REAL(p, N_T_average + 2)*abs(k2 - k1);
}
real liquid_Density(Tracked_Particle *p, real T)
{
	Material *cond_mix = P_MATERIAL(p);
	double Den1, Den2, Den_Droplet;

	if (cond_mix->component[0]->name[0] == 'e')
	{
		
		Den1 = 986.5303588 - 0.6014966034*T - 0.2754046133e-3*T*T;//ace
		Den2 = 1053.6-0.925*T;//eth
		Message("Den1=%f Den2=%f", Den1, Den2);
	}
	else {
		
		Den2 = 986.5303588 - 0.6014966034*T - 0.2754046133e-3*T*T;//ace
		Den1 = 1053.6 - 0.925*T;//eth
	}
	return Den_Droplet = P_USER_REAL(p, N_T_average + 1)*Den1 + P_USER_REAL(p, N_T_average + 2)*Den2;
	Message("Droplet.density=%f", Den_Droplet);
}
real liquid_viscosity(Tracked_Particle *p, real T)
{
	Material *cond_mix = P_MATERIAL(p);
	double k1, k2;

	if (cond_mix->component[0]->name[0] == 'e')
	{
		
		k1 = 0.3183313525e-2-0.1629735179e-4*T+0.223333e-7*T*T; //ace
		k2 = pow(10.0, 686.64/T-5.282); //eth
	}
	else{
		k2 = 0.3183313525e-2 - 0.1629735179e-4*T + 0.223333e-7*T*T;
		k1 = pow(10.0, 686.64 / T - 5.282);
		
		
	}
	return P_USER_REAL(p, N_x_0)*k1 + P_USER_REAL(p, N_x_0 + 1)*k2;
}
real Liquid_diffusivity(Tracked_Particle *p, real T)
{
	Material *cond_mix = P_MATERIAL(p);
	double k1, k2;

//	if (cond_mix->component[0]->name[0] == 'e')
//	{
		
		k1 = 7.4e-12*sqrt(1.0*58.08)*T / (0.3183313525e-2 - 0.1629735179e-4*T + 0.223333e-7*T*T)*1.0e-3*pow(P_USER_REAL(p, N_LJ_l + 1) / 1.18, -0.6 * 3); //Dae 
		k2 = 7.4e-12*sqrt(1.5*46.069)*T / pow(10.0, 686.64 / T - 5.282)*1.e-3*pow(P_USER_REAL(p, N_LJ_l) / 1.18, -0.6 * 3); //Dea 
		//Message("k1=%e k2=%e", k1, k2);
//	}
//	else{
//		k2 = 7.4e-12*sqrt(1.5*46.07)*T / pow(10.0, 686.64 / T - 5.282)*1.e-3*pow(P_USER_REAL(p, N_LJ_l) / 1.18, -0.6 * 3); //ethanol
//		k1 = 7.4e-12*sqrt(1.0*58.08)*T / (0.3183313525e-2 - 0.1629735179e-4*T + 0.223333e-7*T*T)*1.e-3*pow(P_USER_REAL(p, N_LJ_l + 1) / 1.18, -0.6 * 3); //aceton
//	}
	Message("K1=%e K2=%e\n", k1, k2);
		
		return P_USER_REAL(p, N_mdot + 3)*k1 + P_USER_REAL(p, N_mdot + 4 )*k2;
	//Message("K_tot=%f\n", K_tot); swapped
	
}
int molar_surface(Tracked_Particle *p)
{
    real M_tot = 0.e-15;
	
	M_tot = P_USER_REAL(p, 2*N_INT + 1) / 58.08 + P_USER_REAL(p, 3*N_INT + 2) / 46.07;
	
    P_USER_REAL(p, N_mdot + 3) = P_USER_REAL(p, 2*N_INT + 1) / 58.08 / M_tot; //ace
    P_USER_REAL(p, N_mdot + 4) = P_USER_REAL(p, 3*N_INT + 1 + 1) / 46.07 / M_tot; //eth
    return 0;
}

DEFINE_DPM_SCALAR_UPDATE(Diesel_droplet, cell, thread, initialize, p)
{
	int nc = TP_N_COMPONENTS(p);
   	Material *cond_mix = P_MATERIAL(p);

	real Tp = P_T(p);
    real y1, y2;
    real M_tot = 0.e-15;
    int ret;
	
	

	cphase_state_t *c = &(p->cphase);

	clock_t t;
		t = clock();
	//Message("I'm in scalar update\n");
    y1 = p->injection->component.massfr[0];
    y2 = p->injection->component.massfr[1];
	Message("Mfr: y1:%f\t y2:%f\n",y1,y2);
	if (initialize)
	{
		for (int i = 0; i < N_INT + 1; i++) P_USER_REAL(p, i) = Tp;
		for (int i = 0; i < N_INT + 1; i++) P_USER_REAL(p, N_INT + 1 + i) = y1;
		for (int i = 0; i < N_INT + 1; i++) P_USER_REAL(p, 2*N_INT + 2 + i) = y2;
		P_USER_REAL(p, N_T_average) = Tp; //T_av
		P_USER_REAL(p, N_T_average + 1) = y1; //T_av
		P_USER_REAL(p, N_T_average + 2) = y2; //T_av
		P_USER_REAL(p, N_LJ_l) = 4.6; // acetone
		P_USER_REAL(p, N_LJ_l + 1) = 4.53; // ethanol
		P_USER_REAL(p, N_LJ_l + 2) = 560.2; // acetone
		P_USER_REAL(p, N_LJ_l + 3) = 362.6; // ethanol
		P_USER_REAL(p, N_mdot + 5) = 0.0;
		P_USER_REAL(p, N_mdot + 6) = 0.0;
   
	    ret = molar_surface(p);
        ret = gamma_activity(p,Tp);

		//Message("Temperature initiated\n");
		if(cond_mix->component[1]->name[0]=='e') Message("Test names: %c\t %c\n", cond_mix->component[0]->name[0], cond_mix->component[1]->name[0]);
	}
	else
	{
	//	P_USER_REAL(p, 4 * nc + 7 + N_INT + 4) = ((real)t) / CLOCKS_PER_SEC - P_USER_REAL(p, 4 * nc + 7 + N_INT + 3);
		p->state.temp = P_USER_REAL(p, N_T_average); //!return
        //Message("I'm in scalar update\n");
	}

}
DEFINE_DPM_HEAT_MASS(multivap_conv_diffusion_ae, p, Cp, hgas, hvap, cvap_surf, Z, dydt, dzdt)
{
	int ns;
	real Bt = 0.;
	real tot_vap_rate = 0.;
	real vap_rate = 0.;
	int nc = TP_N_COMPONENTS(p);
	Thread *t0 = P_CELL_THREAD(p);
	Material *gas_mix = THREAD_MATERIAL(DPM_THREAD(t0, p));
	Material *cond_mix = P_MATERIAL(p);
	cphase_state_t *c = &(p->cphase);  /* continuous phase struct, caching variables of the cell */
	real Tp;
	real mp = P_MASS(p);
	/* molecular weight in bulk gas */
	real molwt_bulk = 0.;

	//real Dp = DPM_DIAM_FROM_VOL(mp / P_RHO(p));
	real Dp = P_DIAM(p);
	real Ap = DPM_AREA(Dp);
	real Re,Re_liq, Pr, Nu, Sh, Sc, Scl, Nu_star, Sh_Star;
	real Ys, Y_inf, Ys_tot, rho_gas_s, xs_tot, xsM_tot;
	real BM, BT, BT_i, FBM, FBT;
	real x_surf, P_sat, Visc_l, k_l, C_pl;
	real kgas;
	real mugas;
	real T_ref;
	int i, ret, j;
	real c_p_vap;
	Material * cond_c = MIXTURE_COMPONENT(cond_mix, 0);
	real D;
	real L_eff;
	real convective_heating_rate = 0.;
	real T_av;
	real rel_vel, Pe, k_eff;
	real h0, zeta, kappa;
	real lambda[N_Lambda];
	real series[N_Lambda];
	real I_n, b_n;
	real T_eff, dif, coef, phi, dh_dt, h, factor;
	int gas_index;
	real khi_kl = 1.0;
	real khi_kl2 = 1.0;
	real Lat_heat[2];
	//real Velocity = 12.81 - 0.316*(P_TIME(p) / 1000.0 + 0.7884); //100% Acetone
	//real Velocity = 12.30 - 0.344*(P_TIME(p) / 1000.0 + 0.815); //100% Ethanol
	//real Velocity = 12.75 - 0.370*(P_TIME(p) /1000.0 + 0.79); //25% Ethanol 
	//real Velocity = 12.71 - 0.448*(P_TIME(p) / 1000.0 + 0.7952); //50% Ethanol 
	real Velocity = 12.28 - 0.306*(P_TIME(p) / 1000.0 + 0.8233); //75% Ethanol

	//for mass fraction distribution:
	real h0_m, lambda_m0, D_l, D_eff, b_0, s_0, I_0;
	real lambda_m[N_Lambda_y];
	real series_m[N_Lambda_y];
	real temporary; 

	//  real molwt_vap, molwt_vap;


	Message("\ntime(ms)= %f\n", P_TIME(p));
	Tp = P_USER_REAL(p, N_INT);
	T_ref = (c->temp + 2.0*Tp) / 3.0;
	ret = molar_surface(p);

	kgas = Wilke_vap_thermcond(p, T_ref);
	mugas = Wilke_vap_visc(p, T_ref);
	//Message("kgas= %e mugas= %e\n", kgas , mugas);
	
	//Re = p->Re; //?
	//Re = c->rho*Dp*Velocity / c->mu;
	//Message("muc=%e Rhoc=%f\n", c->mu, c->rho);
	Pr = c->sHeat * mugas / kgas; //?
	
	/* Heat transfer coefficient */

	/* molecular weight of gas species */
	/* molwt_bulk is first the reciproke of the real molecular weight to avoid additional divisions */
	mixture_species_loop_i(gas_mix, ns)
	{
		molwt_bulk += c->yi[ns] / solver_par.molWeight[ns]; /* c->yi[ns] contains mass fraction of fluid cell */
	}
	molwt_bulk = MAX(molwt_bulk, DPM_SMALL);
	molwt_bulk = 1. / molwt_bulk;

	/* when not using Runge Kutta solver, increase limiting time for the next integration step */
	if (!p->in_rk)
		p->limiting_time = P_DT(p)*1.01;

	T_av = P_USER_REAL(p, N_T_average);

	ret = gamma_activity(p, T_av);


	/* loop over all components of multi component particle */
	Ys = 0.e-15;
	Y_inf = 0.e-15;
	Ys_tot = 0.e-15;
	xs_tot = 0.e-15;
	xsM_tot = 0.e-15;
	
	Lat_heat[0] = 489.0e3*pow((508.1 - Tp) / (508.1 - MATERIAL_PROP(cond_mix->component[0], PROP_boil_temp)), 0.38);//acetone
	Lat_heat[1] = 120.91e3*pow(516.2 - Tp, 0.38);//ethanol
	
	//Lat_heat[1] = 120.91e3*pow(516.2 - Tp, 0.38);
	//Lat_heat[0] = 394164;
	//Lat_heat[1] = 855237;
	
	//Message("L.Heat Acetone=%f Ethanol=%f\n",Lat_heat[0],Lat_heat[1]);
	//Message("Actual L_Eff= %f\n", L_eff);
	//molar fractions on the surface of the droplet
	for (ns = 0; ns < nc; ns++)
	{
		/* gas species index of vaporization */
		gas_index = TP_COMPONENT_INDEX_I(p, ns);
		if (gas_index >= 0)
		{
			x_surf = P_USER_REAL(p, N_mdot + 3 + ns)*P_USER_REAL(p, N_x_0 + N_components + ns)*exp(MATERIAL_PROP(MIXTURE_COMPONENT(gas_mix, gas_index), PROP_mwi)*Lat_heat[ns] / R_UNIVERSAL*(1.0 / MATERIAL_PROP(cond_mix->component[ns], PROP_boil_temp) - 1.0 / Tp));
			//			x_surf = P_USER_REAL(p, N_mdot + 3 + ns)*P_USER_REAL(p, N_x_0 + N_components + ns)*exp(MATERIAL_PROP(MIXTURE_COMPONENT(gas_mix, gas_index), PROP_mwi)*MATERIAL_PROP(cond_mix->component[ns], PROP_latent_heat) / R_UNIVERSAL*(1.0 / MATERIAL_PROP(cond_mix->component[ns], PROP_boil_temp) - 1.0 / Tp));
			P_USER_REAL(p, N_x_0 + ns) = x_surf;
			xs_tot += x_surf*solver_par.molWeight[gas_index];
			xsM_tot += x_surf;
		}
	}
	xs_tot += (1.0 - xsM_tot)*28.967; //air

	//mass fractions on the surface of the droplet
	for (ns = 0; ns < nc; ns++)
	{
		/* gas species index of vaporization */
		gas_index = TP_COMPONENT_INDEX_I(p, ns);
		if (gas_index >= 0)
		{
			//Message("M_i = %f\n", solver_par.molWeight[gas_index]);
			//Ys = x_surf / ( x_surf + (1.-x_surf) * molwt_bulk / solver_par.molWeight[gas_index]);
			Ys = P_USER_REAL(p, N_x_0 + ns)* solver_par.molWeight[gas_index] / xs_tot;
			//Message("\nGas index = %d\n", gas_index);
			//Message("xs_tot = %f\n", xs_tot);
			//Message("mol weight = %f\n", solver_par.molWeight[gas_index]);
			Y_inf += c->yi[gas_index];
			Ys_tot += Ys;
			P_USER_REAL(p, N_x_0 + 2 * N_components + ns) = Ys;
			P_USER_REAL(p, N_mdot + 5 + ns) = Ys;
		}
	}
	//	if (Ys_tot<1.e-6) L_eff = 489.0*pow((508.1 - Tp) / (508.1 - 329.22), 0.38)*1000.0;
	//	else 	L_eff = L_eff / Ys_tot;
	for (ns = 0; ns < nc; ns++) P_USER_REAL(p, N_x_0 + 2 * N_components + ns) = P_USER_REAL(p, N_x_0 + 2 * N_components + ns) / Ys_tot;

	L_eff = LatentHeat(p, T_av);
	c_p_vap = vapour_heatcap(p, T_ref); //vapour ?
	D = vapour_diffusivity(p, T_ref, c->pressure);
	rho_gas_s = Vapour_Density(p, T_ref, c->pressure);//!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Message("vapour density %f\n",rho_gas_s);
	Sc = mugas / (rho_gas_s * D);//?
	//Message("Vap diffusivity = %e\n c_p_vap = %f\nL_eff = %f, Sc= %f\n", , c_p_vap, L_eff,Sc);
	//Message("\nVapour density= %f", rho_gas_s);
	Message("Actual L_Eff= %f\n",L_eff);
	//if (Ys_tot<1.e-6) L_eff = 489.0*pow((508.1 - Tp) / (508.1 - 329.22), 0.38)*1000.0;
	//else 	L_eff = L_eff / Ys_tot;
	//Re = rho_gas_s*Dp*Velocity / c->mu;
	Re = rho_gas_s*Dp*Velocity / mugas;
	
	//Message(" muc=%e Rhoc=%f\n", c->mu, c->rho);


	BM = (Ys_tot - Y_inf) / (1.0 - Ys_tot);
	//BM = (Ys_tot) / (1.0 - Ys_tot);
	FBM = 1.0;
	if (BM > 1.e-12) 
	FBM = pow(1.0 + BM, 0.7)*log(1.0 + BM) / BM;;
	Sh_Star = 2.0 + (pow(1.0 + Re*Sc, 1.0 / 3.0)*MAX(1.0, pow(Re, 0.077)) - 1.0) / FBM;
	Sh = log(1.0 + BM)*Sh_Star;
	//Sh = log(1.0 + BM)*(2.0 + 0.6*sqrt(Re)*pow(Sc, 1.0 / 3.0));
	tot_vap_rate = Ap * D * rho_gas_s * Sh / Dp; //?
	P_USER_REAL(p, N_mdot) = tot_vap_rate;
	P_USER_REAL(p, N_mdot + 1) = P_USER_REAL(p, N_x_0 + 2 * N_components)*tot_vap_rate;
	P_USER_REAL(p, N_mdot + 2) = P_USER_REAL(p, N_x_0 + 2 * N_components + 1)*tot_vap_rate;


	BT = BM;
	BT_i = BT;
	dif = 1.0;
	coef = c_p_vap * rho_gas_s * D / kgas * Sh_Star;//?
	//Message(" Vap.Cp= %f\nVap.Diff= %e\t Vap.Ther.cond= %f\n",c_p_vap,D,kgas);

	phi = 0.e-15;
	while (dif > ACCURACY)
	{
		FBT = 1.0;
		if (BT > 1.0e-12) 
		FBT = pow(1.0 + BT, 0.7)*log(1.0 + BT) / BT;
		Nu_star = 2.0 + (pow(1.0 + Re*Pr, 1.0 / 3.0)*MAX(1.0, pow(Re, 0.0777)) - 1.0) / FBT;
		phi = coef / Nu_star;
		BT = pow(1.0 + BM, phi) - 1.0;
		dif = abs(BT - BT_i);
		BT_i = BT;
	}

	Nu = Nu_star;
	if (BT > 1.0e-12) 
	Nu = log(1.0 + BT) * Nu_star / BT;

	//Message("\nBM= %f BT= %f Phi= %f\nNu= %f Sh= %f Nu_star= %f Sh star= %f\nRe= %f Pr= %f coef= %f\n",BM,BT,phi,Nu,Sh,Nu_star,Sh_Star,Re,Pr,coef);
	//Message("\nBM= %f FBM= %f BT= %f FBT=%f\nNu= %f Nu_star= %f Sh= %f Sh star= %f\nRe= %f Pr= %f coef= %f\n", BM, FBM, BT, FBT, Nu, Nu_star, Sh, Sh_Star, Re, Pr, coef);
	Visc_l = liquid_viscosity(p, T_av); //?
	k_l = liquid_thermal_conductivity(p, T_av);
	Re_liq = P_RHO(p)*Velocity * Dp / Visc_l;
	//k_l *= 1000.0;
	C_pl = p->Cp;

	//Message("\nLiq.viscosity=%e\t liq.therm.cond=%f\t CONST=%f liq.Cp=%f\n",Visc_l,k_l,khi_kl,C_pl);

	//rel_vel = sqrt((c->V[0] - P_VEL(p)[0])*(c->V[0] - P_VEL(p)[0]) + (c->V[1] - P_VEL(p)[1])*(c->V[1] - P_VEL(p)[1]) + (c->V[2] - P_VEL(p)[2])*(c->V[2] - P_VEL(p)[2]));
	rel_vel = Velocity;
	Pe = 12.69 / 16.0*P_RHO(p)*0.5*Dp* C_pl / k_l*rel_vel* mugas / Visc_l*pow(Re, 1.0 / 3.0) / (1.0 + BM);
	
	if (abs(Pe) > 1.e-12)
		
	khi_kl = 1.86 + 0.86*tanh(2.225*log10(Pe / 30.0));
	k_eff = khi_kl*k_l;

	//T_eff = c->temp - tot_vap_rate*L_eff / PI / Dp / Nu / kgas;
	T_eff = c->temp - (tot_vap_rate*L_eff) / (PI * Dp * Nu * kgas);
	h0 = kgas*Nu*0.5 / k_eff - 1.0;
	zeta = (h0 + 1.0)*T_eff;//?
	kappa = k_eff / (C_pl*P_RHO(p)*0.25*Dp*Dp);
	
	//Message("\nKappa= %f rel_vel = %f Visc_l= %e k_eff=%f\nPe = %f Cpl= %f dropD=%f zeta = %f\n",kappa, rel_vel, Visc_l, k_eff, Pe, C_pl, P_RHO(p), zeta);

	//temperature distribution
	ret = Lambda(h0, lambda);
	for (i = 0; i < N_Lambda; i++) series[i] = 0.e-15;

	for (i = 0; i < N_Lambda; i++)
	{
		b_n = 0.5*(1.0 + h0 / (h0*h0 + lambda[i] * lambda[i]));
		I_n = 0.e-15;
		I_n = P_USER_REAL(p, N_INT)*sin(lambda[i]);
		for (j = 1; j < N_INT; j += 2) I_n += 4.0 * P_USER_REAL(p, j)*(((double)j)*Delta_R)*sin(lambda[i] * ((double)j)*Delta_R);
		for (j = 2; j < N_INT; j += 2) I_n += 2.0 * P_USER_REAL(p, j)*(((double)j)*Delta_R)*sin(lambda[i] * ((double)j)*Delta_R);
		I_n = I_n*Delta_R / 3.0;
		series[i] = (I_n - sin(lambda[i]) / lambda[i] / lambda[i] * zeta)*exp(0.0 - kappa*lambda[i] * lambda[i] * P_DT(p)) / b_n;
		
		//Message("bn= %f In=%f",b_n,I_n);
	}

	for (j = 0; j < N_INT + 1; j++) P_USER_REAL(p, j) = T_eff;
	for (i = 0; i < N_Lambda; i++)
	{
		P_USER_REAL(p, 0) += series[i] * lambda[i];
		for (j = 1; j < N_INT + 1; j++) P_USER_REAL(p, j) += series[i] * sin(lambda[i] * ((double)j)*Delta_R) / (((double)j)*Delta_R);
	}


	Tp = P_USER_REAL(p, N_INT);

	T_av = Tp;
	for (j = 1; j < N_INT; j += 2) T_av += 4.0 * P_USER_REAL(p, j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
	for (j = 2; j < N_INT; j += 2) T_av += 2.0 * P_USER_REAL(p, j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
	T_av = T_av*Delta_R;
	//Message("T_eff= %f Tav= %f\n",T_eff,T_av);

	//mass fraction distribution acetone
	D_l = Liquid_diffusivity(p, T_av);
	Scl = Visc_l / P_RHO(p) / D_l;
	//Sc = Visc_l / P_RHO(p) * D_l;
	if (abs(Pe) > 1.e-12) 
	khi_kl2 = 1.86 + 0.86*tanh(2.225*log10(Re_liq * Scl / 30.0));
	D_eff = khi_kl2*D_l;
	h0_m = -tot_vap_rate* 0.5 / (D_eff * PI * Dp * P_RHO(p)) - 1.0;
	
	//Message("\nLiq.Diff= %e Sc = %e Scl=%e\nD_eff = %e Const=%f\n Rp=%e alpha = %e\n", D_l, Sc,Scl, D_eff, khi_kl, Dp/2.0, tot_vap_rate/PI/P_RHO(p)/Dp/Dp);
	//Message("\nKappa= %f rel_vel = %f\t k_eff=%f\nPe = %f h0 = %f\t zeta = %f\n", kappa, rel_vel, k_eff, Pe, h0, zeta);
	//Message("\nLiq.Diff= %e D_eff = %e\n", D_l, D_eff);
	//Message("h0 = %f\t h0_m = %f l0 = %f\n", h0, h0_m, lambda_m0);
	//Message("h0_m = %f\t d_eff=%f\n",h0_m, D_eff);
	//Message("Pdens=%e Dp=%e\n", P_RHO(p), Dp);
	
	Message("Pdens=%e", P_RHO(p));
	lambda_m0 = -1.0;
	lambda_m0 = Lambda_m0(h0_m);

	for (i = 0; i < N_Lambda_y; i++) lambda_m[i] = -1.0;
	ret = Lambda(h0_m, lambda_m);

	for (i = 0; i < N_Lambda_y; i++) series_m[i] = 0.e-15;

	
	for (i = 0; i < N_Lambda; i++)
	
	//Message("l0 = %f\n", lambda_m0);
	//Message("lambda = %f\t", lambda_m[i]);
	//Message("\n");
	//for (j = 0; j < N_INT + 1; j++) Message("Y1[%d] = %e\n", j, P_USER_REAL(p, N_INT + 1 + j));


	b_0 = -0.5*(1.0 + h0_m / (h0_m*h0_m - lambda_m0 * lambda_m0));
	I_0 = 0.e-15;
	I_0 = P_USER_REAL(p, N_INT + 1 + N_INT)*sinh(lambda_m0);
	for (j = 1; j < N_INT; j += 2) I_0 += 4.0 * P_USER_REAL(p, N_INT + 1 + j)*(((double)j)*Delta_R)*sinh(lambda_m0 * ((double)j)*Delta_R);
	for (j = 2; j < N_INT; j += 2) I_0 += 2.0 * P_USER_REAL(p, N_INT + 1 + j)*(((double)j)*Delta_R)*sinh(lambda_m0 * ((double)j)*Delta_R);
	I_0 = I_0*Delta_R / 3.0;
	s_0 = (I_0 + P_USER_REAL(p, N_x_0 + 2 * N_components)*sinh(lambda_m0) / lambda_m0 / lambda_m0 * (1.0 + h0_m))*exp(4.0*D_eff*lambda_m0 * lambda_m0 / Dp / Dp * P_DT(p)) / b_0;

	//Message("b_0= %e I_0= %e s_0= %e", b_0, I_0, s_0);

	for (i = 0; i < N_Lambda_y; i++)
	{
		b_n = 0.5*(1.0 + h0_m / (h0_m*h0_m + lambda_m[i] * lambda_m[i]));
		I_n = 0.e-15;
		I_n = P_USER_REAL(p, N_INT + 1 + N_INT)*sin(lambda_m[i]);
		for (j = 1; j < N_INT; j += 2) I_n += 4.0 * P_USER_REAL(p, N_INT + 1 + j)*(((double)j)*Delta_R)*sin(lambda_m[i] * ((double)j)*Delta_R);
		for (j = 2; j < N_INT; j += 2) I_n += 2.0 * P_USER_REAL(p, N_INT + 1 + j)*(((double)j)*Delta_R)*sin(lambda_m[i] * ((double)j)*Delta_R);
		I_n = I_n*Delta_R / 3.0;

		//Message("b_n= %e I_n= %e s_n= %e", b_n, I_n, series_m[i]);

		series_m[i] = (I_n - P_USER_REAL(p, N_x_0 + 2 * N_components)*sin(lambda_m[i]) / lambda_m[i] / lambda_m[i] * (1.0 + h0_m))*exp(0.0-4.0*D_eff*lambda_m[i] * lambda_m[i] / Dp / Dp * P_DT(p)) / b_n;
	}
	//Message("b_n= %e I_n= %e s_n= %e", b_n, I_n, series_m[i]);

	for (j = 0; j < N_INT + 1; j++) P_USER_REAL(p, N_INT + 1 + j) = P_USER_REAL(p, N_x_0 + 2 * N_components);
	//Message("\nY_1[0] = %f ", P_USER_REAL(p, N_INT + 1));
	P_USER_REAL(p, N_INT + 1) += s_0 * lambda_m0;
	//Message("Y_1[0] = %f ", P_USER_REAL(p, N_INT + 1));
	for (j = 1; j < N_INT + 1; j++) P_USER_REAL(p, N_INT + 1 + j) += s_0 * sinh(lambda_m0 * ((double)j)*Delta_R) / (((double)j)*Delta_R);
	for (i = 0; i < N_Lambda_y; i++)
	{
		P_USER_REAL(p, N_INT + 1) += series_m[i] * lambda_m[i];
		//Message("Y_1[0] = %f ", P_USER_REAL(p, N_INT + 1));
		for (j = 1; j < N_INT + 1; j++) P_USER_REAL(p, N_INT + 1 + j) += series_m[i] * sin(lambda_m[i] * ((double)j)*Delta_R) / (((double)j)*Delta_R);
	}
	//Message("Y_1[0] = %f \n", P_USER_REAL(p, N_INT + 1));

	//for (j = 0; j < N_INT + 1; j++) P_USER_REAL(p, N_INT + 1 + j) = 1.0; //acetone
	
	//Message("b_0= %e I_0= %e s_0= %e", b_0, I_0, s_0);
	//Message("b_n= %e I_n= %e s_n= %e", b_n, I_n, series_m[i]);
	//for (i = 0; i < N_Lambda_y; i++) series_m[i] = 0.e-15;
	//Message("l0 = %f\n", lambda_m0);
	//for (i = 0; i < N_Lambda; i++) Message("lambda = %f\t", lambda_m[i]);
	//Message("\n");


	//Mass fraction ethanol
	for (i = 0; i < N_Lambda_y; i++) series_m[i] = 0.e-15;

	b_0 = -0.5*(1.0 + h0_m / (h0_m* h0_m - lambda_m0 * lambda_m0));
	I_0 = 0.e-15;
	I_0 = P_USER_REAL(p, 2 * N_INT + 2 + N_INT)*sinh(lambda_m0);
	for (j = 1; j < N_INT; j += 2) I_0 += 4.0 * P_USER_REAL(p, 2 * N_INT + 2 + j)*(((double)j)*Delta_R)*sinh(lambda_m0 * ((double)j)*Delta_R);
	for (j = 2; j < N_INT; j += 2) I_0 += 2.0 * P_USER_REAL(p, 2 * N_INT + 2 + j)*(((double)j)*Delta_R)*sinh(lambda_m0 * ((double)j)*Delta_R);
	I_0 = I_0*Delta_R / 3.0;
	s_0 = (I_0 + P_USER_REAL(p, N_x_0 + 2 * N_components + 1)*sinh(lambda_m0) / lambda_m0 / lambda_m0 * (1.0 + h0_m))*exp(4.0*D_eff*lambda_m0 * lambda_m0 / Dp / Dp * P_DT(p)) / b_0;
	
	//Message("\nb_0= %e I_0= %e s_0= %e\n", b_0, I_0, s_0);

	for (i = 0; i < N_Lambda_y; i++)
	{
		b_n = 0.5*(1.0 + h0_m / (h0_m* h0_m + lambda_m[i] * lambda_m[i]));
		I_n = 0.e-15;
		I_n = P_USER_REAL(p, 3*N_INT + 2)*sin(lambda_m[i]);
		for (j = 1; j < N_INT; j += 2) I_n += 4.0 * P_USER_REAL(p, 2 * N_INT + 2 + j)*(((double)j)*Delta_R)*sin(lambda_m[i] * ((double)j)*Delta_R);
		for (j = 2; j < N_INT; j += 2) I_n += 2.0 * P_USER_REAL(p, 2 * N_INT + 2 + j)*(((double)j)*Delta_R)*sin(lambda_m[i] * ((double)j)*Delta_R);
		I_n = I_n*Delta_R / 3.0;
		series_m[i] = (I_n - P_USER_REAL(p, N_x_0 + 2 * N_components + 1)*sin(lambda_m[i]) / lambda_m[i] / lambda_m[i] * (1.0 + h0_m))*exp(0.0-4.0*D_eff*lambda_m[i] * lambda_m[i] / Dp / Dp * P_DT(p)) / b_n;
	}

	for (j = 0; j < N_INT + 1; j++) P_USER_REAL(p, 2 * N_INT + 2 + j) = P_USER_REAL(p, N_x_0 + 2 * N_components + 1);


	P_USER_REAL(p, 2 * N_INT + 2) += s_0 * lambda_m0;
	for (j = 1; j < N_INT + 1; j++) P_USER_REAL(p, 2 * N_INT + 2 + j) += s_0 * sinh(lambda_m0 * ((double)j)*Delta_R) / (((double)j)*Delta_R);
	for (i = 0; i < N_Lambda_y; i++)
	{
		P_USER_REAL(p, 2 * N_INT + 2) += series_m[i] * lambda_m[i];
		for (j = 1; j < N_INT + 1; j++) P_USER_REAL(p, 2 * N_INT + 2 + j) += series_m[i] * sin(lambda_m[i] * ((double)j)*Delta_R) / (((double)j)*Delta_R);
	}

	//for (j = 0; j < N_INT + 1; j++) P_USER_REAL(p, 2 * N_INT + 2 + j) = 0.0;//ethanol

	//Normalisation y0 + y1 = 1
	//	for (j = 0; j < N_INT + 1; j++) if (P_USER_REAL(p, N_INT + 1 + j) >1.0) P_USER_REAL(p, N_INT + 1 + j) = 1.0;
	//	for (j = 0; j < N_INT + 1; j++) if (P_USER_REAL(p, N_INT + 1 + j) < 1.e-5) P_USER_REAL(p, N_INT + 1 + j) = 0.0;

	//	for (j = 0; j < N_INT + 1; j++) if (P_USER_REAL(p, 2 * N_INT + 2 + j) >1.0) P_USER_REAL(p, 2 * N_INT + 2 + j) = 1.0;
	//	for (j = 0; j < N_INT + 1; j++) if (P_USER_REAL(p, 2 * N_INT + 2 + j) < 1.e-5) P_USER_REAL(p, 2 * N_INT + 2 + j) = 0.0;

	//	for (j = 0; j < N_INT + 1; j++)
	//{
	//temporary = P_USER_REAL(p, N_INT + 1 + j)+ P_USER_REAL(p, 2 * N_INT + 2 + j);
	//P_USER_REAL(p, N_INT + 1 + j) = P_USER_REAL(p, N_INT + 1 + j) / temporary;
	//P_USER_REAL(p, 2 * N_INT + 2 + j) = P_USER_REAL(p, 2 * N_INT + 2 + j) / temporary;
	//}
	
	

	P_USER_REAL(p, N_T_average + 1) = P_USER_REAL(p, N_INT + 1 + N_INT);
	for (j = 1; j < N_INT; j += 2) P_USER_REAL(p, N_T_average + 1) += 4.0 * P_USER_REAL(p, N_INT + 1 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
	for (j = 2; j < N_INT; j += 2) P_USER_REAL(p, N_T_average + 1) += 2.0 * P_USER_REAL(p, N_INT + 1 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
	P_USER_REAL(p, N_T_average + 1) = P_USER_REAL(p, N_T_average + 1)*Delta_R;
	
	//Message("Y1_av= %f", P_USER_REAL(p, N_T_average + 1));


	P_USER_REAL(p, N_T_average + 2) = P_USER_REAL(p, 2 * N_INT + 2 + N_INT);
	for (j = 1; j < N_INT; j += 2) P_USER_REAL(p, N_T_average + 2) += 4.0 * P_USER_REAL(p, 2*N_INT + 2 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
	for (j = 2; j < N_INT; j += 2) P_USER_REAL(p, N_T_average + 2) += 2.0 * P_USER_REAL(p, 2*N_INT + 2 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
	P_USER_REAL(p, N_T_average + 2) = P_USER_REAL(p, N_T_average + 2)*Delta_R;
	P_USER_REAL(p, N_T_average + 2) = 1.0 - P_USER_REAL(p, N_T_average + 1);
	for (j = 0; j < N_INT + 1; j++) P_USER_REAL(p, N_INT + N_INT + 2 + j) = 1.0 - P_USER_REAL(p, N_INT + 1 + j);
	
	//Message("\nb_0= %e I_0= %e s_0= %e\n", b_0, I_0, s_0);
	//Message("b_n= %e I_n= %e s_n= %e\n", b_n, I_n, series_m[i]);
	
	Message("Y1_av= %f", P_USER_REAL(p, N_T_average + 1));
	Message("Y2_av= %f\n", P_USER_REAL(p, N_T_average + 2));
	
	for (ns = 0; ns < nc; ns++)
	{
		/* gas species index of vaporization */
		gas_index = TP_COMPONENT_INDEX_I(p, ns);
		if (gas_index >= 0)
		{
			vap_rate = P_USER_REAL(p, N_x_0 + 2 * N_components + ns) * tot_vap_rate; //
//			vap_rate = tot_vap_rate; //!!
//			if (Ys_tot>1.e-6) vap_rate = vap_rate*P_USER_REAL(p, nc + ns) / Ys_tot; //!!

			if (!p->in_rk)
				if (ABS(vap_rate)>0.)
					p->limiting_time = MIN(p->limiting_time, dpm_par.fractional_change_factor_mass*P_MASS(p) / vap_rate*TP_COMPONENT_I(p, ns));

			P_USER_REAL(p, N_mdot + 1 + ns) = vap_rate;
			//dydt[1 + ns] -= vap_rate;
			dydt[1 + ns] -= 0.0;

			{
				int source_index = injection_par.yi2s[gas_index];
				if (source_index >= 0)
				{
					dzdt->species[source_index] += vap_rate;
					//p->source.mtc[source_index] = c->rho * Ap * Sh_Star * D / Dp;
					p->source.mtc[source_index] = c->rho * PI * Dp * Sh_Star * D;//?
				}
			}
		}
	}


	//p->source.htc = Nu*kgas*P_DIAM(p);
	p->source.htc = 0.e-15;
	dh_dt = Nu * kgas * Ap / Dp * (c->temp - T_av);//?
	//Message("dh_dt = %f\n", Nu * kgas * Ap / Dp * (c->temp - T_av));

	h = Nu * kgas / Dp;
	convective_heating_rate = h * Ap / (mp * p->Cp);

	//dydt[0] += dh_dt / (mp * p->Cp);
	//dzdt->energy -= dh_dt;
	/* limit for higher heating rate */
	if (!p->in_rk)
	{
		if (ABS(convective_heating_rate)>DPM_SMALL)
		{
			factor = dpm_par.fractional_change_factor_heat;
			if (ABS(c->temp - Tp)>Tp)
				factor = dpm_par.fractional_change_factor_heat*Tp / (c->temp - Tp);
			p->limiting_time = MIN(p->limiting_time, factor / ABS(convective_heating_rate));
		}
	}
	P_USER_REAL(p, N_T_average) = T_av;
	p->state.temp = T_av;
	dydt[0] = 0.e-15;
}
/*-------------------Gas mixture properties-------------------------*/
DEFINE_PROPERTY(ethanol_vapour_ThermCond, c, t)
{
	real temp = C_T(c, t);
	return 1.8037e-3+7.419e-6*temp+1.1536e-7*temp*temp;//Sazhin et al (2010) IJHMT
}
DEFINE_PROPERTY(acetone_vapour_ThermCond, c, t)
{
	real temp = C_T(c, t);
	return 0.01143468*pow(temp / 273.15, 2.0 - 0.342843365e-6*temp*temp);//Sazhin et al (2010) IJHMT
}
DEFINE_PROPERTY(acetone_vapour_viscosity, c, t)
{
	real temp = C_T(c, t);
	return -0.1579884447e-5+0.3090958433e-7*temp-0.449329e-11*temp*temp;//Sazhin et al (2010) IJHMT
}
DEFINE_PROPERTY(ethanol_vapour_viscosity, c, t)
{
	real temp = C_T(c, t);
	return -0.19757e-6 + 0.29211e-7*temp;//Sazhin et al (2010) IJHMT
}
DEFINE_PROPERTY(Wilke_thermcond, c, t)
{
	real sum = 0.;
	real phi[9];
	real visc[3], molwt[3], tc[3],x[3];
	real Phi, M_tot;
	int i,j;
	Material *spi;
	Material *spj;
	Material *gas_mix = THREAD_MATERIAL(t);
	Property *prop;
	real result = 0.e-15;
	
    real T = C_T(c,t);

	M_tot = 0.e-15;

    visc[0] = -0.19757e-6+0.29211e-7*T; //ethanol
    visc[1] = -0.1579884447e-5 + 0.3090958433e-7*T-0.449329e-11*T*T; //acetone
    visc[2] = 1.8e-5;

    molwt[0] = 46.07; //ethanol
    molwt[1] = 58.08; //acetone
    molwt[2] = 28.967;

    tc[0] = 1.8037e-3+7.419e-6*T+1.1536e-7*T*T; //vap thermcond ethanol
    tc[1] = 0.01143468*pow(T/273.15,2.0-0.342843365e-6*T*T); //vap thermcond acetone
    tc[2] = 0.0242;

	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		//prop = (MATERIAL_PROPERTY(spi));
		//visc[i] = generic_property(c, t, prop, PROP_mu, C_T(c, t));
		//tc[i] = generic_property(c, t, prop, PROP_ktc, C_T(c, t));
		//molwt[i] = MATERIAL_PROP(MIXTURE_COMPONENT(gas_mix, i), PROP_mwi);
		M_tot += C_YI(c, t, i) / molwt[i];
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			if (i == j) phi[Ind(i, j)] = 1.0;
			else phi[Ind(i, j)] = (1.0+sqrt(visc[i]/visc[j]*sqrt(molwt[j]/molwt[i]))) / sqrt(8.0) / (sqrt(1.0 + molwt[i] / molwt[j]));
		}
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		x[i] = C_YI(c, t, i) / molwt[i] / M_tot;
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		Phi = 0.e-15;
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			Phi += x[j] * phi[Ind(i, j)];
		}
		result += x[i] * tc[i] / Phi;
	}


	return result;
}
DEFINE_PROPERTY(Wilke_visc, c, t)
{
	real sum = 0.;
	real phi[9];
	real visc[3], molwt[3], x[3];
	real Phi, M_tot;
	int i, j;
	Material *spi;
	Material *spj;
	Material *gas_mix = THREAD_MATERIAL(t);
	Property *prop;
	real result = 0.e-15;
    real T = C_T(c,t);
	
	M_tot = 0.e-15;

    visc[0] = -0.19757e-6+0.29211e-7*T; //ethanol
    visc[1] = -0.1579884447e-5 + 0.3090958e-7*T-0.449329e-11*T*T; //acetone
    visc[2] = 1.8e-5;

    molwt[0] = 46.07; //ethanol
    molwt[1] = 58.08; //acetone
    molwt[2] = 28.967;
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		//prop = (MATERIAL_PROPERTY(spi));
		//visc[i] = generic_property(c, t, prop, PROP_mu, C_T(c, t));
		//molwt[i] = MATERIAL_PROP(MIXTURE_COMPONENT(gas_mix, i), PROP_mwi);
		M_tot += C_YI(c, t, i) / molwt[i];
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			if (i == j) phi[Ind(i, j)] = 1.0;
			else phi[Ind(i, j)] = (1.0 + sqrt(visc[i] / visc[j] * sqrt(molwt[j] / molwt[i]))) / sqrt(8.0) / (sqrt(1.0 + molwt[i] / molwt[j]));
		}
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		x[i] = C_YI(c, t, i) / molwt[i] / M_tot;
	}
	mixture_species_loop(THREAD_MATERIAL(t), spi, i)
	{
		Phi = 0.e-15;
		mixture_species_loop(THREAD_MATERIAL(t), spj, j)
		{
			Phi += x[j] * phi[Ind(i, j)];
		}
		result += x[i] * visc[i] / Phi;
	}
	
	return result;
	
}
/*-----------------Droplet liquid properties-------------------------*/
DEFINE_DPM_PROPERTY(Diesel_liquid_density, c, t, p, T)
{
	return 744.11 - 0.771*(P_T(p) - 300.0);
	//or set constant density
}
DEFINE_DPM_PROPERTY(Diesel_liquid_specific_heat, c, t, p, T)
{
	return (2.18 + 0.0041*(P_T(p) - 300.0))*1000.0;
}
DEFINE_DPM_PROPERTY(Acetone_Latent_Heat, c, t, p, T)
{
	return 489.0 * pow((508.1 - T) / (508.1 - 329.22), 0.38)*1000.0;
}
DEFINE_DPM_PROPERTY(Ethanol_Latent_Heat, c, t, p, T)
{
	return 120.91e3*pow(516.2-T,0.38);
}
DEFINE_DPM_PROPERTY(Droplet_Specific_Heat, c, t, p, T)
{
	return 2165.234225 - 2.963*P_T(p) + 0.01*P_T(p)*P_T(p);
}
DEFINE_DPM_PROPERTY(Droplet_Density, c, t, p, T)
{
	//return  P_USER_REAL(p, N_T_average + 1)*Den1 + P_USER_REAL(p, N_T_average + 2)*Den2;
	
}
DEFINE_DPM_TIMESTEP(Constant_dt, p, dt)
{
	return 1.e-6;
}
DEFINE_DPM_PROPERTY(Mixture_liquid_specific_heat, c, t, p, T)
{
	int nc = TP_N_COMPONENTS(p);
	int ns;
	real result = 0.e-15;
	Material *cond_mix = P_MATERIAL(p);
	//Message("Mass fractions: 0: %f\t; 1: %f\n", P_USER_REAL(p, N_T_average + 1), P_USER_REAL(p, N_T_average + 2));
	//Message("C_pl fractions: 0: %f\t; 1: %f\n", MATERIAL_PROPERTY(cond_mix->component[0], PROP_Cp), MATERIAL_PROP(cond_mix->component[1], PROP_Cp));
    //Message("C_pl fractions: 0: %f\n", cond_mix->component[0]->p->udf_prop .fcn());
	//Message("N_components %d\n", nc);
	/*for (ns = 0; ns < nc; ns++)
	{
		result += P_USER_REAL(p, N_T_average + 1 + ns)*MATERIAL_PROP(cond_mix->component[ns], PROP_Cp);//Y_i*C_li
	}
	return result;*/
    return P_USER_REAL(p, N_T_average + 1)*(2165.234225-2.963*T+0.01*T*T)+P_USER_REAL(p, N_T_average + 2)*(15039.0-130.53*T+0.4143*T*T-0.39583e-3*T*T*T);
}
DEFINE_DPM_PROPERTY(Mixture_liquid_Density, c, t, p, T)
{
	int nc = TP_N_COMPONENTS(p);
	int ns;
	real result = 0.e-15;
	Material *cond_mix = P_MATERIAL(p);
	
	return P_USER_REAL(p, N_T_average + 1)*(986.5303588 - 0.6014966034*T - 0.2754046133e-3*T*T) + P_USER_REAL(p, N_T_average + 2)*(1053.6 - 0.925*T);
}

DEFINE_DPM_PROPERTY(Acetone_liquid_specific_heat, c, t, p, T)
{
	return 2165.234225-2.963*T+0.01*T*T;
}
DEFINE_DPM_PROPERTY(Ethanol_liquid_specific_heat, c, t, p, T)
{
	return 15039.0-130.53*T+0.4143*T*T-0.39583e-3*T*T*T;
}
DEFINE_DPM_PROPERTY(Acetone_liquid_density, c, t, p, T)
{
	return 986.5303588-0.6014966034*T-0.2754046133e-3*T*T;

}
DEFINE_DPM_PROPERTY(Ethanol_liquid_density, c, t, p, T)
{
	return 1053.6-0.925*T;
}