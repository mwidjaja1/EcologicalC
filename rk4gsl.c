//
//  main.c
//  Rk4gsl
//
//  Created by Matthew Widjaja on 11/19/13.
//  Copyright (c) 2013 Matthew Widjaja. All rights reserved.
//

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <complex.h>

int func (double t, const double y[], double f[], void *params)
{
	const int maxNode = 8;
	
	double newAlpha[maxNode][maxNode] = {
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	};
	
	double constR[maxNode] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	
	double steadyValP[maxNode] = {
		38.8089080995067, 43.5963099130963, 39.4611963199623,
		27.6303237065429, 26.8372058084276, 28.0030748854752,
		34.7552952146415, 47.0617038749445
	};
	
	double oldAlpha[maxNode][maxNode] = {
		{-2.3305, -4.7368, -0.8303, -1.4475, 0, 0, 0, -0.291868997002398},
		{-0.2565, -4.7009, -0.522, -0.3013, 0, 0, 0, -0.647544780938589},
		{-0.1342, -1.1164, -1.7607, -0.5695, 0, 0, 0, -0.523688100581974},
		{-0.2881, -0.3794, -0.2665, -2.968, 0, 0, 0, -1.01620423202169},
		{0, 0, 0, -0.478736026755033, -3.1091, -0.0561, -1.0181, -0.8356},
		{0, 0, 0, -1.08243797034983, -0.0955, -2.1481, -1.1953, -0.347},
		{0, 0, 0, -0.242649494404053, -0.8154, -0.4252, -2.7833, -0.576},
		{0, 0, 0, -1.12489771091476, -0.3553, -0.0653, -0.3707, -2.3896}
	};
	
	for (int i1=0; i1<maxNode; i1++) {
		for (int i2=0; i2<maxNode; i2++) {
			newAlpha[i1][i2] = -oldAlpha[i1][i2] / steadyValP[i1];
			constR[i1] = constR[i1] + (newAlpha[i1][i2] * steadyValP[i2]);
		}
	}
	
	
	for (int i1=0; i1<maxNode; i1++) {
		f[i1] = constR[i1] * y[i1];
		for (int i2=0; i2<8; i2++) {
			f[i1] = f[i1] - (newAlpha[i1][i2] * f[i1] * f[i2]);
		}
	}
	
	return GSL_SUCCESS;
}


int jac ()
{
	return GSL_SUCCESS;
}


int main (void)
{
	const int maxNode = 8;
	double maxTime = 200.0;
	
	// Declares the System
	gsl_odeiv2_system sys = {func, jac, maxNode};
	
	
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new
		(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
	
	double t = 0.0, t1 = maxTime;
	double y[maxNode] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	
	for (int i=1; i<=maxTime; i++)
	{
      double ti = i * t1 / maxTime;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
		
      if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		
		printf ("Time: %f\n",t);
		printf ("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
				  y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
		
	}
	
	gsl_odeiv2_driver_free (d);
	return 0;
}