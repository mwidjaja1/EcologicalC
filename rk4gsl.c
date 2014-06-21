/* --- MATTHEW WIDJAJA: LOTKA-VOLTERA MODEL ------------------------------------
 ** Purpose: This C code solves a series of Lotka-Volterra ODEs across a period
 **			 of time to solve how various organisms interact in a Food Web.
 **
 ** Credit:
 **	Matthew Widjaja			Umm, well... I wrote this. This is awkward.
 **	Russ Manson				For assistance in using C & the GSL Library.
 **	Jason Shulman			For assistance in the original MATLAB algorithm.
 **	GNU Scientific Library	Provided the ODE Library (libgsl0-dev) used.
 **	Michael Widjaja			Reminded me why I should have never grown up.
 **
 ** Compiling Rule: (Note, replace /opt/local/lib w. /usr/local/lib in Ubuntu)
 **		gcc -L/opt/local/lib rk4gsl.c -lgsl -lgslcblas -lm -std=c99
 **
 ** This code consists of four functions:
 **	Func A	Equations	This GSL-mandated function creates the ODE functions
 **						with the user-specified constants for solving.
 **	Func B	Jacobson	This GSL-mandated function is a placeholder because
 **						the RK45 solver being utilized does not require this.
 **	Func C	ODE Pass	This function takes the variables of Func D & passes
 **						it to Func A and works with Func A to solve the ODEs.
 **						This function then passes the results back to Func D.
 **	Func D	Main		This function is the primary program of this script.
 **						The user selects how they want to test the ODEs here &
 **						the proper modifications for those queries are made.
 **
 ** The globally defined variables & parameters (used in Func A, C, & D) are:
 **	maxNode	const int	States the max amount of nodes present
 **	fixQty	int			States the amount of fixed nodes
 **	fixEqu	int			Makes a [maxNode] sized array to identify fixed nodes
 **	svFile	File		Stores all Y & T values if requested in Func D.
 **
 ** These variables might have to be modified prior to running the model:
 **	maxNode		In Global	States the max amount of nodes present
 **	steadyValP	In Func A	Final steady-state values expected
 **	oldAlpha	In Func A	The old Alpha Values
 **	maxTime		In Func D	The max time the ODE Solver should calculate for
 **
 ** Unless otherwise specified in Func D, all output data is stored in:
 **	outFile		In Func D	This file stores the last set of each run's data.
 **	odeResult	In Func C	This [maxNode+1] sized array stores the last set
 **							of each run's data.
 **	sumResult	In Func D	This [maxRun][maxNode+1] sized array stores the
 **							the last set of data from every run. The script
 **							automatically builds this, if deemed necessary.
 **	svFile		In Func C	This file stores all Y & T values from each run.
 **							This is optional: User selects this in Func D.2.
 ** ------------------------------------------------------------------------- */


// Import of Function Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


// Global Functions used for Functions A, C, & D
#define maxNode 8
int fixQty;
int fixEqu[maxNode];
FILE* svFile;


/* --- FUNCTION A: EQUATIONS ---------------------------------------------------
 ** This function generates the Equations that will be solved by the ODE Solver.
 ** 
 ** Many of the parameteters & variables used in Func D were declared earlier
 **	as global commands. These are: maxNode, fixQty, and fixEqu.
 **
 ** The parameters & variables for this include:
 **	newAlpha	double	A [maxNode] x [maxNode] sized array to hold the new
 **						alpha values that this function will later produce.
 **	constR		double	A [maxNode] sized array to hold the constant values
 **						that this function will later create.
 **	steadyValP	double	A [maxNode] sized array with the steady state values
 **						each node should reach at the end of this model.
 **	oldAlpha	double	A [maxNode] x [maxNode] sized array w. old alpha
 **						values. This is the unified birth & death rates.
 **
 ** The Algorithm used for this function is:
 **	1. Declare all of the parameters & variables aforementioned.
 **	2. Designate a forLoop to create the newAlpha values & constR values.
 **	3. Specify the Differential Equations using the values obtained in A.2.
 **	3A. Checks if a node was fixed (as determined by Func D) and if so, does
 **		not make its function.
 **	3B. If 3A was not true, then the Node's Rates & Function are created.
 ** ------------------------------------------------------------------------- */
int func (double t, const double y[], double f[], void* params)
{
	
	// --- Step A.1: Declare the Parameters & Variables ------------------------
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
	
	// --- Step A.2: Generate the newAlpha & constR values ---------------------
	for (int i1=0; i1<maxNode; i1++) {
		for (int i2=0; i2<maxNode; i2++) {
			newAlpha[i1][i2] = -oldAlpha[i1][i2] / steadyValP[i1];
			constR[i1] = constR[i1] + (newAlpha[i1][i2] * steadyValP[i2]);
		}
	}
	
	// --- Step A.3: Create the ODE's Di5fferential Equations ------------------
	int fixNode=0;
	for (int i1=0; i1<maxNode; i1++) {
		
		// --- Step A.3A: Checks if Node was Fixed -----------------------------
		for (int i2=0; i2<fixQty; i2++) {
			if (fixEqu[i2] == i1) {
				f[i1] = 0;
				fixNode = 1;
			}
		}
		
		// --- Step A.3B: Creates Rates & Functions ----------------------------
		if (fixNode == 0) {
			f[i1] = constR[i1] * y[i1];
			for (int i3=0; i3<maxNode; i3++) {
				f[i1] = f[i1] - (newAlpha[i1][i3] * y[i1] * y[i3]);
			}
		}
		fixNode = 0;
	}
	
	return GSL_SUCCESS;
}	// Function A Terminates Here



/* --- FUNCTION B: JACOBSON ----------------------------------------------------
 ** This function would normally declare the Jacobson Matrix. The RK45 solver
 ** used does not require this. Even though this function isn't actually needed,
 ** GSL still requires such a Jacobson function to be declared with no value.
 ** ------------------------------------------------------------------------- */
int jac ()
{
	return GSL_SUCCESS;
}	// Function B Terminates Here



/* --- FUNCTION C: ODE_PASS ----------------------------------------------------
 ** This function solves the ODEs with an RK45 method. 
 **
 ** These parameters & variables needed for this function are declared in
 ** Func D (the Main Function) and are linked to this one. They include:
 **	t			double		The initial time.
 **	t1			double		The final time (when the model should stop).
 **	maxTime		double		This should be the same value as 't1'.
 **	y			double		The [maxNode] sized array which holds final results.
 **	intRun		int			Tracks how often the ODE Solver was ran.
 **
 ** If svInt = 2 (as documented in Func D), the following variables are used:
 **	svBuffer	char[25]	Creates the file name for each output file made.
 **	svFile		FILE		Declared globally, this is the file that is used.
 **
 ** The Algorithm used for this function is:
 **	1. We first declare the ODE system & the ODE driver which GSL requires.
 **	2. If svInt = 2 (ie. the user wanted all data posted on the terminal &
 **    saved to a file), each run's output file is created here.
 **	3. A forLoop is used to state how frequently the results should be saved.
 **		3A. The stepsize is calculated
 **		3B. The ODE solver is ran
 **		3C. (If svInt = 1 or 2) The results are displayed to the terminal.
 **		3D. (If svInt = 2) The results are saved to the output file
 **	4. The final results are saved to an array which is linked to Func D.
 **	   The ODE solver is then shut down & (if svInt = 2) svFile is saved.
 ** ------------------------------------------------------------------------- */
double odePass (double t, double t1, double maxTime, double y[], int intRun,
					 double odeResult[], int svInt) {
	
	// --- Step C.1: Declare ODE System & Driver -------------------------------
	gsl_odeiv2_system sys = {func, jac, maxNode};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new
		(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
	
	// --- Step C.2: Declare Output File Parameters ----------------------------
	if (svInt==2) {		// Runs if svInt==2 (it means the user wanted all data)
		char svBuffer[25];
		sprintf(svBuffer, "output_%d.dat", intRun);
		svFile = fopen(svBuffer, "w");
	
		if(NULL == svFile)
			printf("\n-No allFile at intrun=%i. Pick bigger stepSize-", intRun);
	}
	
	// --- Step C.3: States how often should Results be Saved ------------------
	for (int i=1; i<=maxTime; i++) {
      
		// --- Step C.3A: Calculate Stepsize -----------------------------------
		double ti = i * t1 / maxTime;
		
		// --- Step C.3B: ODE Solver -------------------------------------------
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}
		
		if (svInt==1 || svInt == 2) {
		// --- Step C.3C: If svInt = 1 or 2, Displays Results to Terminal ------
		printf ("\nTime: %f\n",t);
		printf ("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
				  y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
		}
		
		if (svInt == 2) {
		// --- Step C.3D: If svInt = 2, Saves Results to an Output File --------
		fprintf (svFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
					t, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
		}
		
	}
	
	
	// --- Step C.4: Stores & Saves Final Results ------------------------------
	// Populates the Answer Array
	odeResult[0] = intRun;
	for(int i=1; i<maxNode+1; i++)
		odeResult[i] = y[i-1];
	return *odeResult;				// Sends the Answer Array to Func D
	
	gsl_odeiv2_driver_free (d);	// Terminates the ODE Driver/Solver
	
	if (svInt==2)	// Closes the allFile if user wanted all data
		fclose(svFile);
	
}	// Function C Terminates Here



/* --- FUNCTION D: MAIN_PROGRAM ------------------------------------------------
 ** This program is the 'main' program of this script. This controls how the
 ** equations are modified before processing by the ODE solver.
 **
 ** These variables, defined here in Func D, will primarily be used in Func C:
 **	t, t1, maxTime, maxNode, y, intRun, svInt.
 **
 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
 **
 ** DO NOT MODIFY. These parameters & variables are strictly used in Func D
 ** for all Func D queries:
 **	intMethod	int		Saves the method the user selected in Step D.2.
 **	yInitial	double	A [maxNode] sized array with Y's initial values.
 **	y			double	A [maxNode] sized array with Y's initial values,
 **						after the user modifies its values.
 **
 ** The svInt interger in Step D.2 decides how data is saved. Large svInt =
 ** longer processing times & a greater likelihood of memory crashes.
 **	svInt = 0	Only the last set of data from each run is saved into 
 **				outFile, odeResult, & (if applicable) sumResult.
 **	svInt = 1	svInt=0 applies. Also, all Y & T data from each run will be
 **				displayed on the Terminal Prompt.
 **	svInt = 2	svInt=0 and svInt=1 applies. Also, all Y & T data from each
 **				run will be saved to an output file, as documented in Func C.
 **
 ** The Algorithm used for this function is:
 **	1. We declare the parameters & variables which Functions D (and C) need.
 **	2. We create the output array & file which will store the final set of
 **	   data from every run. (svInt = 0, 1, or 2)
 **	3. The user gets to select how they want data to be saved. (via svInt)
 **	4. The user gets to select a query -- this are options on how data can be
 **	   obtained. (Documentation for these queries are noted before each one).
 ** ------------------------------------------------------------------------- */
int main (int argc, char **argv) {

	// --- Step D.1: Declare Parameters & Variables ----------------------------
	int intMethod = 0;	// Asks user for the query to be used.
	int svInt = 0;			// Asks if user wants all Y & T values saved.
	int intRun = 0;		// States the Run Number the program is on
	
	double t = 0.0;								// Initial Time for Model
	double maxTime = 200, t1 = maxTime;			// Max Time for the Model
	
	double yInit[maxNode] =						// Initial Conditions for Y
	{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double y[maxNode];							// Creates Array for Modified Y
	
	for (int i=0; i<maxNode; i++)				// Sets Initial Y = Modified Y
		y[i] = yInit[i];
	
	
	// --- Step D.2: Creation of Output Array & File ---------------------------
	double odeResult[maxNode+1];
	FILE *outFile = fopen("results.dat", "w");
	fprintf(outFile, "Run Number\t Y0\t Y1\t Y2\t Y3\t Y4\t Y5\t Y6\t Y7\n");

	
	// --- Step D.3: Asks how Data should be Saved -----------------------------
	printf("--Rk4gsl is now running--\n");
	printf("\n--Output Options--\n");	// Selects if all data should be saved
		printf("  0. (Default) Only save the last set of data of each run.\n");
		printf("  1. Do #0 & output all Y & T results to the Terminal\n");
		printf("  2. Do #0 & #1 & save #1 results to a file\n");
		printf("\nSelect an Output Option (0-2): ");
		scanf("%d",&svInt);
	
	
	// --- Step D.4: Selection of Query ----------------------------------------
	printf("\n--Query Options--\n");	// Selects Query
		printf("  0. Manual\n  1. Automated SKO\n");
		printf("  2. Automated DKO\n  3. Automated TKO\n");
		printf("  4. Single Node Sweep\n  5. Double Node Sweep\n");
		printf("\nSelect a Query (0-5): ");
		scanf("%d",&intMethod);
	printf("\n--Query %d is Running--\n",intMethod);	// Confirms the Query
	
	
	/* --- QUERY D.4.0: Manual ODE ---------------------------------------------
	 ** This query lets the user run the ODE Function & fix specified nodes.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 ** These Func D defined variables will be used: intRun, outFile, odeResult
	 **
	 ** The following parameters & variables will be modified in Query 0:
	 **	fixVal		User Set	States the value each node should be fixed with.
	 **
	 ** The following parameter & variables are unique to Query 0:
	 **	fixNode		User Set	States which nodes should be fixed.
	 ** --------------------------------------------------------------------- */
	if (intMethod == 0) {
		// Asks for how many nodes to fix
		int fixNode = 0;
		double fixVal;
		printf("How many nodes should be fixed? = "); scanf("%d",&fixQty);
		
		// Asks for which nodes to fix & to which values
		if (fixQty >= 1) {
			printf("\nThe range of nodes is 0 to %d",maxNode-1);
			for (int i=0; i<fixQty; i++) {
				printf("\nKnockout Node Number= "); scanf("%d",&fixNode);
				fixEqu[i] = fixNode;
				printf("Knockout with the Value = "); scanf("%lf",&fixVal);
				y[fixNode] = fixVal;
				intRun++;
			}
		}
		
		// Runs Function C -- the ODE Solver
		odePass(t, t1, maxTime, y, intRun, odeResult, svInt);
		
		// Saves results to an output file & saves the file
		fprintf (outFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
					odeResult[0], odeResult[1], odeResult[2], odeResult[3],
					odeResult[4], odeResult[5], odeResult[6], odeResult[7],
					odeResult[8]);
		fclose(outFile);

	}


	/* --- QUERY D.4.1: Automated SKO ------------------------------------------
	 ** This query lets the user do an automated Single Knock Out routine.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 ** These Func D defined variables will be used: intRun, outFile, odeResult
	 **
	 ** The following parameters & variables will be modified in Query 1:
	 **	fixVal		Script Set		States the value a node will be fixed with.
	 ** --------------------------------------------------------------------- */
	if (intMethod == 1) {
		fixQty = 1;		// SKO -- This fixes 1 node at a time
		double sumResult[maxNode][maxNode+1];		// Creates results array
		
		for (int i1=0; i1<maxNode; i1++) {
			fixEqu[0] = i1; y[i1] = 0;		// Fixes each Node Automatically
			odePass(t, t1, maxTime, y, intRun, odeResult, svInt);	// Func C
			y[i1] = yInit[i1];				// Resets the Initial Conditions

			// Saves Results to Output Array
			for (int i2=0; i2<maxNode+1; i2++)
				sumResult[intRun][i2] = odeResult[i2];
			
			// Saves Results to Output File
			fprintf (outFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
						sumResult[intRun][0], sumResult[intRun][1],
						sumResult[intRun][2], sumResult[intRun][3],
						sumResult[intRun][4], sumResult[intRun][5],
						sumResult[intRun][6], sumResult[intRun][7],
						sumResult[intRun][8]);
			
			intRun++;	// Adds 1 to the intRun Counter

		}
	}
	
	
	/* --- QUERY D.4.2: Automated DKO ------------------------------------------
	 ** This query lets the user do an automated Double Knock Out routine.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 ** These Func D defined variables will be used: intRun, outFile, odeResult
	 **
	 ** The following parameters & variables will be modified in Query 2:
	 **	fixVal	Script Set		States the value a node will be fixed with.
	 **
	 ** The following parameters & variables are unique to Query 2:
	 **	maxRun	const int		Specifies the maximum amount of runs possible.
	 ** --------------------------------------------------------------------- */
	if (intMethod == 2) {
		fixQty = 2;								// DKO -- Fixes 2 nodes at once
		
		const int maxRun = pow(maxNode,2);		// Determines max amount of runs
		double sumResult[maxRun][maxNode+1];	// Creates the Solutions Array
		
		// Fixes Node 1
		for (int i1=0; i1<maxNode; i1++) {
			fixEqu[0] = i1; y[i1] = 0;
			
			// Fixes Node 2
			for (int i2=0; i2<maxNode; i2++) {
				// If both nodes are same, we skip it.
				if (i1==i2) {i2++; continue;}
				
				// Otherwise, we solve both nodes
				fixEqu[1] = i2; y[i2] = 0;	// Fixes the Second Node
				odePass(t, t1, maxTime, y, intRun, odeResult, svInt); // Func C
				y[i2] = yInit[i2];			// Resets Node 2 to original value
				
				// Saves Results to Output Array
				for (int i3=0; i3<maxNode+1; i3++)
					sumResult[intRun][i3] = odeResult[i3];
				
				// Saves Results to Output File
				fprintf (outFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
							sumResult[intRun][0], sumResult[intRun][1],
							sumResult[intRun][2], sumResult[intRun][3],
							sumResult[intRun][4], sumResult[intRun][5],
							sumResult[intRun][6], sumResult[intRun][7],
							sumResult[intRun][8]);
				
				intRun++;	// Adds 1 to intRun Counter
			}
			y[i1] = yInit[i1];	// Resets Node 1 to its Original Value
		}
	}
	
	
	/* --- QUERY D.4.2: Automated TKO ------------------------------------------
	 ** This query does an automated Triple Knock Out Routine.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 ** These Func D defined variables will be used: intRun, outFile, odeResult
	 **
	 ** The following parameters & variables will be modified in Query 2:
	 **	fixVal	Script Set		States the value a node will be fixed with.
	 **
	 ** The following parameters & variables are unique to Query 2:
	 **	maxRun	const int		Specifies the maximum amount of runs possible.
	 ** --------------------------------------------------------------------- */
	if (intMethod == 3) {
		fixQty = 3;		// TKO -- This fixes 3 nodes at a time
		
		const int maxRun = pow(maxNode,fixQty);	// Determines max amount of runs
		double sumResult[maxRun][maxNode+1];	// Creates the Solutions Array
		
		// Fixes Node One
		for (int i1=0; i1<maxNode; i1++) {
			fixEqu[0] = i1; y[i1] = 0;
			
			// Fixes Node Two
			for (int i2=0; i2<maxNode; i2++) {
				if (i1==i2) {i2++; continue;}	// Skip Node 2 if Node2=Node1
				fixEqu[1] = i2; y[i2] = 0;		// Otherwise, fix the 2nd Node
				
				// Fixes Node Three
				for (int i3=0; i3<maxNode; i3++) {
					// If Nodes 3 & 1 OR 3 & 2 are same, we skip said Node 3
					if (i3==i1 || i3==i2) {i3++; continue;}
					
					// Otherwise, we solve for all three nodes
					fixEqu[2] = i3; y[i3] = 0;	// Fixes the First Node
					odePass(t, t1, maxTime, y, intRun, odeResult, svInt);
					y[i3] = yInit[i3];		// Resets Node 3 to original value
				
					// Saves Results to Output Array
					for (int i4=0; i4<maxNode+1; i4++)
						sumResult[intRun][i4] = odeResult[i4];
				
					// Saves Results to Output File
					fprintf (outFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
							 sumResult[intRun][0], sumResult[intRun][1],
							 sumResult[intRun][2], sumResult[intRun][3],
							 sumResult[intRun][4], sumResult[intRun][5],
							 sumResult[intRun][6], sumResult[intRun][7],
							 sumResult[intRun][8]);
				
				intRun++;	// Adds 1 to intRun Counter
				}
				y[i2] = yInit[i2];	// Resets Node 2 to its original value
			}
			y[i1] = yInit[i1];		// Resets Node 1 to its original value
		}
	}
	

	/* --- QUERY D.4.4: Automated One-Node Sweep -------------------------------
	 ** This query lets the user do an automated sweep of one node as its value
	 ** changes from 0 to a max value in specified stepSize increments.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 ** These Func D defined variables will be used: intRun, outFile, sumResult,
	 **												 odeResult
	 **
	 ** The following parameters & variables will be modified in Query 4:
	 **	fixVal		Script Set	States the value a node will be fixed with.
	 **
	 ** The following parameters & variables are unique for Query 4:
	 **	fixNode		User Set		States which nodes should be fixed (in 4B).
	 **	maxVal		User Set		States the max value the fixed node will
	 **								reach (in 4B).
	 **	stepSize	User Set		States how often Func C will be ran (in 4C).
	 **
	 ** The Algorithm used for Query 4 is:
	 **	4A. Creates Parameters & Variables for the code.
	 **	4B. Asks the users on which node to fix and to what maximum value.
	 **	4C. Allocates memory for the solutions array.
	 **	4D. A For-Loop is designated to sweep node 1 from 0 to its max value.
	 **	4E. A 2nd For-Loop is designated to retrieve steady-state data from
	 **		Func C and save it to the output file (results.dat).
	 ** --------------------------------------------------------------------- */
	if (intMethod == 4) {
		// --- Step D.4.4A: Declaring Parameters & Variables -------------------
		fixQty = 1;		// SKO -- This fixes 1 node at a time
		int fixNode = 0;
		double maxVal, stepSize;
		
		// --- Step D.4.4B: Asks for which nodes to fix to which max values ----
		printf("Fix Node # (Note: The 1st node is node 0) = ");
			scanf("%d",&fixNode);
			fixEqu[0] = fixNode;
		printf("State the Fixed Node's Maximum Value = "); scanf("%lf",&maxVal);
				
		// --- Step D.4.4C: Allocates Memory for Solutions ---------------------
		printf("Select the Stepsize to be used = "); scanf("%lf",&stepSize);
		const int maxRun = maxVal/stepSize;
		double sumResult[maxRun][maxNode+1];
		
		// --- Step D.4.4D: Sweeps Node 1 from 0 -> its Max Value --------------
		for (double i=0; i<=maxVal; i=i+stepSize) {
			y[fixNode] = i;		// Fixes the Node to the current 'Sweep' value
			odePass(t, t1, maxTime, y, intRun, odeResult, svInt);	// Func C
			printf("Node %d = %f\n",fixNode,i);
			
			// --- Step D.4.4E: Retrieves Data & Outputs to File ---------------
			for (int i2=0; i2<maxNode; i2++)
				sumResult[intRun][i2] = odeResult[i2+1];
			fprintf (outFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
						sumResult[intRun][0], sumResult[intRun][1],
						sumResult[intRun][2], sumResult[intRun][3],
						sumResult[intRun][4], sumResult[intRun][5],
						sumResult[intRun][6], sumResult[intRun][7],
						sumResult[intRun][8]);

			intRun++;	// Notes that a new iteration will begin
		}
		fclose(outFile);	// Saves the Output File

	}
	
	
	/* --- QUERY D.4.5: Automated Two-Node Sweep -------------------------------
	 ** This query lets the user do an automated sweep of two nodes as each of
	 ** the nodes have its values approach from 0 to its specified max value
	 ** (which could be different for both) in a specified stepSize increment.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 ** These Func D defined variables will be used: intRun, outFile, sumResult,
	 **												 odeResult
	 **
	 ** The following parameters & variables will be modified in Query 5:
	 **	fixVal		Script Set	States the value a node will be fixed with.
	 **
	 ** The following parameters & variables are unique for Query 5:
	 **	fixNode		User Set	States which nodes should be fixed.
	 **	tempVal		User Set	States the max value the fixed node will reach.
	 **	maxVal		Script Set	A [fixQty] sized array storing the max values
	 **							all nodes will reach.
	 **	stepSize	User Set	States how often Func C should be ran.
	 **
	 ** The Algorithm used for Query 5 is:
	 **	5A. Creates Parameters & Variables for the code.
	 **	5B. Asks the users on which nodes to fix and to what maximum value.
	 **	5C. Allocates memory for the solutions array. NOTE: If an output file
	 **		 will be saved for each Y/T value from each run, a large stepSize 
	 **		 should be selected or the code might crash from excessive memory.
	 **	5D. A For-Loop is made to sweep node 1 from 0 to its max value.
	 **	5E. A 2nd For-Loop is made to sweep node 2 from 0 to its max value.
	 **	5F. A 3rd For-Loop is made to retrieve steady-state data from
	 **		 Func C and save it to the output file (results.dat).
	 ** ---------------------------------------------------------------------- */
	if (intMethod == 5) {
		// --- Step D.4.5A: Declaring Parameters & Variables -------------------
		fixQty = 2;		// DKO -- This fixes 2 nodes at a time
		int fixNode = 0;
		double tempVal, maxVal[2], stepSize;
		
		// --- Step D.4.5B: Asks for which nodes to fix to which max values ----
		for (int i=0; i<2; i++) {	// This loop automates the query
			printf("Fix Node # (Note: The 1st node is node 0) = ");
				scanf("%d",&fixNode);
				fixEqu[i] = fixNode;
			printf("Select the Max Value of the Fixed Node = ");
				scanf("%lf",&tempVal);
				maxVal[i] = tempVal;
		}
		
		// --- Step D.4.5C: Allocates Memory for Solutions ---------------------
		if (svInt == 2)	// Warns the user of potential crashes if applicable
			printf("\nBecause svInt=2, a small Stepsize (<1) may crash!\n");
		 
		printf("Select the Stepsize to be used = "); scanf("%lf",&stepSize);
		const int maxRun = 1.5 * (maxVal[0]/stepSize) * (maxVal[1]/stepSize);
		double sumResult[maxRun][maxNode+1];
		
		// --- Step D.4.5D: Sweeps Node 1 from 0 -> its Max Value --------------
		for (double i1=0; i1<=maxVal[0]; i1=i1+stepSize) {
			y[fixEqu[0]] = i1;		// Fixes the First Node Automatically
			printf("Node %d = %f\n",fixEqu[0],i1);
			
			// --- Step D.4.5E: Sweeps Node 2 from 0 -> its Max Value ----------
			for (double i2=0; i2<=maxVal[1]; i2=i2+stepSize) {
				y[fixEqu[1]] = i2;		// Fixes the Second Node Automatically
				odePass(t, t1, maxTime, y, intRun, odeResult, svInt);  // Func C
				
				// --- Step D.4.5F: Retrieves Data & Outputs to File ---------------
				for (int i3=0; i3<=maxNode+2; i3++)
					sumResult[intRun][i3] = odeResult[i3];
				fprintf (outFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
							sumResult[intRun][0], sumResult[intRun][1],
							sumResult[intRun][2], sumResult[intRun][3],
							sumResult[intRun][4], sumResult[intRun][5],
							sumResult[intRun][6], sumResult[intRun][7],
							sumResult[intRun][8]);
				
				intRun++;	// Notes that a new iteration will begin
			}
		}
		fclose(outFile);	// Saves the Output File
	}
	printf("\n--The process has completed successfully--\n");
}	// Function D Terminates Here