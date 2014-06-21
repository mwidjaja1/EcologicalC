//
//  main.c
//  Rk4gsl
//
//  Created by Matthew Widjaja on 11/19/13.
//  Copyright (c) 2013 Matthew Widjaja. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv2.h>
#include "gsl_odeiv_modified.h"



/* --- FUNCTION A: EQUATIONS ---------------------------------------------------
 ** This function generates the Equations that will be solved by the ODE Solver.
 ** 
 ** The parameters & variables for this include:
 ** maxNode		int		The maximum amount of nodes used.
 **							(note, maxNode is also declared seperately in Func D)
 ** newAlpha	double	A [maxNode] x [maxNode] sized array to hold the new
 **							alpha values that this function will later produce.
 ** constR		double	A [maxNode] sized array to hold the constant values
 **							that this function will later create.
 ** steadyValP	double	A [maxNode] sized array with the steady state values
 **							each node should reach at the end of this model.
 ** oldAlpha	double	A [maxNode] x [maxNode] sized array w. old alpha values.
 **
 ** The Algorithm used for this function is:
 ** 1. Declare all of the parameters & variables aforementioned.
 ** 2. Designate a forLoop to create the newAlpha values & constR values.
 ** 3. Specify the Differential Equations using the values obtained in step 2.
 ** ------------------------------------------------------------------------- */
int func (double t, const double y[], double f[], void *params)
{
	
	// --- Step A.1: Declare the Parameters & Variables -------------------------
	const int maxNode = 8;
	double newAlpha[maxNode][maxNode] = {
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	};
	double constR[maxNode] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	double steadyValP[maxNode] = {
		38.8089080995067, 43.5963099130963, 39.4611963199623, 27.6303237065429,
		26.8372058084276, 28.0030748854752, 34.7552952146415, 47.0617038749445
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
	
	// --- Step A.2: Generate the newAlpha & constR values ----------------------
	for (int i1=0; i1<maxNode; i1++) {
		for (int i2=0; i2<maxNode; i2++) {
			newAlpha[i1][i2] = -oldAlpha[i1][i2] / steadyValP[i1];
			constR[i1] = constR[i1] + (newAlpha[i1][i2] * steadyValP[i2]);
		}
	}
	
	// --- Step A.3: Create the ODE's Differential Equations --------------------
	for (int i1=0; i1<maxNode; i1++) {
		f[i1] = constR[i1] * y[i1];
		for (int i2=0; i2<8; i2++) {
			f[i1] = f[i1] - (newAlpha[i1][i2] * y[i1] * y[i2]);
		}
	}
	
	return GSL_SUCCESS;
}



/* --- FUNCTION B: JACOBSON ----------------------------------------------------
 ** This function would normally declare the Jacobson Matrix. The RK45 solver
 ** used does not require this. Even though this function isn't actually needed,
 ** GSL still requires such a Jacobson function to be declared with no value.
 ** ------------------------------------------------------------------------- */
int jac ()
{
	return GSL_SUCCESS;
}



/* --- FUNCTION C: ODE_FUNCTION ------------------------------------------------
 ** This function solves the ODEs with an RK45 method. 
 **
 ** The parameters & variables needed for this function are declared in
 ** Function D (the Main Function) and are linked to this one. They include:
 ** t				double	The initial time
 ** t1			double	The final time (when the model should stop)
 ** maxTime		double	This should be the same value as 't1'
 ** maxNode		int		The maximum amount of nodes used.
 **							(note, maxNode is also declared seperately in Func A)
 ** y				double	The maxNode sized array which holds the final results
 ** OTHERS		Note that the file path for the output file & buffer in Step C.2
 **				must be declared. This file should be a '.dat' file for plotting.
 **
 ** The Algorithm used for this function is:
 ** 1. We first declare the ODE system & the ODE driver which GSL requires.
 ** 2. Declare the Output File & Buffer so that the results can be saved.
 ** 3. A forLoop is declared to state how frequent the ODE solver should be ran.
 ** 3A. The stepsize is calculated
 ** 3B. The ODE solver is ran
 ** 3C. The results are displayed to the terminal & saved to the output file.
 ** ------------------------------------------------------------------------- */
void odeFunc (double t, double t1, double maxTime,
				  const int maxNode, double y[]) {
	
	// --- Step C.1: Declare ODE System & Driver --------------------------------
	gsl_odeiv2_system sys = {func, jac, maxNode};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new
		(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
	
	// --- Step C.2: Declare Output File & Buffer -------------------------------
	char outBuffer[100];
	FILE *out = fopen("/Users/Matthew/Dropbox/Academics/CPLS/Final Project/Rk4gsl/Rk4gsl/output.dat", "w");
	
	for (int i=1; i<=maxTime; i++) {
      // --- Step C.3A: Calculate Stepsize -------------------------------------
		double ti = i * t1 / maxTime;
		
		// --- Step C.3B: ODE Solver ---------------------------------------------
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}
		
		// --- Step C.3C: Display & Save Results ---------------------------------
		printf ("\nTime: %f\n",t);		// Display Results
		printf ("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
				  y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
		
		fputs(outBuffer, out);			// Saves Results
		fprintf (out,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
					t, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
	}
	
	gsl_odeiv2_driver_free (d);
	fclose(out);
}



/* --- FUNCTION D: MAIN_PROGRAM ------------------------------------------------
 ** This program is the 'main' program of this script. This controls how the
 ** equations are modified before processing by the ODE solver.
 **
 ** Many of the parameters & variables declared in Func D are used in Func C, &
 ** thus are explained there. These are: t, t1, maxTime, maxNode, & y.
 **
 ** The parameters & variables that are declared & only used in Func D are:
 ** intMethod		int		Saves the method the user selected in Step D.2.
 ** masterQty		int		The quantity of master nodes.
 ** masterNodes	int		A [masterQty] sized array which states which nodes
 **								are master nodes.
 ** yInitial		double	A [maxNode] sized array with Y's initial values.
 ** y					double	A [maxNode] sized array with Y's initial values,
 **								after the user modifies 
 **
 ** The Algorithm used for this function is:
 ** 1. We first declare the parameters & variables which Functions D & C need.
 ** 2. The user gets to select which query/experiment they'd like to run.
 **	 (Further documentation for thse queries are noted prior to each query)
 ** ------------------------------------------------------------------------- */
int main (int argc, char **argv) {

	// --- Step D.1: Declare Parameters & Variables -----------------------------
	int intMethod = 0;								// States the Method the user picked

	const int maxNode = 8;							// Qty of Nodes
	const int masterQty = 2;						// Qty of Master Nodes
	int masterNodes[masterQty] = {4, 8};		// States nodes are Master Nodes
	
	double t = 0.0;									// Initial Time for Model
	double maxTime = 200, t1 = maxTime;			// Max Time for the Model
	
	double yInitial[maxNode] =						// Initial Conditions for Y
	{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double y[maxNode];								// Creates Array for Modified Y Values
	
	for (int i=0; i<maxNode; i++)					// Sets Initial Y = Modified Y for now
		y[i] = yInitial[i];
	
	// --- Step D.2: Selection of Query -----------------------------------------
	printf("Select a Method: "); scanf("%d",&intMethod);
	printf("IntMethod: %d\n", intMethod);
	
	/* --- QUERY D.A: Manual ODE ------------------------------------------------
	 ** This query permits the user to run the ODE Equations and specify which
	 ** nodes they want knocked out, if any. This is the default option.
	 **
	 ** The parameters & variables declared in this query is:
	 ** fixQty	int		User-specified value of how many nodes should be fixed.
	 ** fixNode	int		User-specified value of which nodes are fixed.
	 ** fixVal	double	User-specified value of a node's initial value.
	 **
	 ** The Algorithm used for this function is:
	 ** A1. The user gets to select how many nodes should be fixed
	 ** A2. An ifStatement is declared if the user said a node should be fixed.
	 **	 In it, a forLoop is declared to have the user specify which nodes
	 **	 should be fixed and with what initial value.
	 ** A3. The odeFunc (FuncC) is ran to produce results.
	 ** ---------------------------------------------------------------------- */
	if (intMethod == 0) {
		// --- Step D.2.A1: Asks for Qty of Nodes to Fix -------------------------
		int fixQty; int fixNode = 0; double fixVal;
		printf("\nHow many nodes should be fixed? = "); scanf("%d",&fixQty);
		
		// --- Step D.2.A2: Asks for which Nodes to Fix & how so -----------------
		if (fixQty >= 1) {
			double fixEqu[fixQty];
			printf("\nThe range of nodes is 0 to %d",maxNode-1);
			for (int i=0; i<fixQty; i++) {
				printf("\nKnockout Node Number= "); scanf("%d",&fixNode);
				fixEqu[i] = fixNode;
				printf("Knockout with the (double) Value = "); scanf("%lf",&fixVal);
				y[fixNode] = fixVal;
			}
		}
		
		// --- Step D.2.A3: Runs ODE Solver --------------------------------------
		odeFunc(t, t1, maxNode, maxTime, y );
	}
	
	
	
}