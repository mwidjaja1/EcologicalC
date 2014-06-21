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

int func (double t, const double y[], double f[], void *params)
{
	f[0] = y[0] * (1 - y[0] - y[1]);
	f[1] = y[1] * (1 - y[1] + y[0]);
	return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
	gsl_matrix_view dfdy_mat
	= gsl_matrix_view_array (dfdy, 0, 0);
	gsl_matrix * m = &dfdy_mat.matrix;
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}

int main (void)
{
	gsl_odeiv2_system sys = {func, jac, 2};
	
	gsl_odeiv2_driver * d =
	gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
											 1e-6, 1e-6, 0.0);
	int i;
	double t = 0.0, t1 = 20.0;
	double y[2] = { 3.0, 1.0 };
	
	for (i = 1; i <= 100; i++)
	{
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
		
      if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		
      printf ("%f\t %f\t %f\n", t, y[0], y[1]);
	}
	
	gsl_odeiv2_driver_free (d);
	return 0;
}