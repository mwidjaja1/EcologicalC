/* --- Matthew Widjaja: Final Project ------------------------------------------
 ** Theoratically, this code solves a series of ODEs across a period of time.
 ** Realistically, this code solves a Lotka-Volterra Style Food Web.
 **
 ** This code uses the GNU Scientific Library for completion.
 ** This code is copyrighted under wImagineering.
 **
 ** This code consists of four functions:
 ** Func A	Equations	This GSL-mandated function creates the ODE functions
 **							with the user-specified constants for solving.
 ** Func B	Jacobson		This GSL-mandated function is a placeholder because the
 **							RK45 solver being utilized does not require this.
 ** Func C	ODE Pass		This function takes the variables of Func D & passes it
 **							over to Func A. This function also saves the results.
 ** Func D	Main			This function is the primary program of this script.
 **							The user selects how they want to test the ODEs here.
 **
 ** The globally defined variables & parameters (used in Func A & D) are:
 ** maxNode	const int	States the max amount of nodes present
 ** fixQty	int			States the amount of fixed nodes
 ** fixEqu	int			Creates a [maxNode] sized array to identify fixed nodes
 **
 */


// Import of Function Libraries
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


// Global Functions used for Functions A, B, C, & D
#define maxNode 8;	// States the max amount of nodes present
#define fixQty 2;			// States the current amount of master nodes
int fixEqu[8];		// Creates the matrix which will identify fixed nodes



/* --- FUNCTION A: EQUATIONS ---------------------------------------------------
 ** This function generates the Equations that will be solved by the ODE Solver.
 ** 
 ** Many of the parameteters & variables used in Func D were declared earlier
 ** as global commands. These are: maxNode, fixQty, and fixEqu.
 **
 ** The parameters & variables for this include:
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
 ** 3A. Checks if a node was fixed (as determined by Func D) and if so, doesn't
 **	  make its function.
 ** 3B. If 3A was not true, then the Node's Rates & Function are created.
 ** ------------------------------------------------------------------------- */
int func (double t, const double y[], double f[], void* params)
{
	
	// --- Step A.1: Declare the Parameters & Variables -------------------------
	double newAlpha[8][8] = {
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
		{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}, {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	};
	double constR[8] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	
	double steadyValP[8] = {
		38.8089080995067, 43.5963099130963, 39.4611963199623, 27.6303237065429,
		26.8372058084276, 28.0030748854752, 34.7552952146415, 47.0617038749445
	};
	double oldAlpha[8][8] = {
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
	for (int i1=0; i1<8; i1++) {
		for (int i2=0; i2<8; i2++) {
			newAlpha[i1][i2] = -oldAlpha[i1][i2] / steadyValP[i1];
			constR[i1] = constR[i1] + (newAlpha[i1][i2] * steadyValP[i2]);
		}
	}
	
	// --- Step A.3: Create the ODE's Differential Equations --------------------
	int fixNode=0;
	for (int i1=0; i1<maxNode; i1++) {
		
		// --- Step A.3A: Checks if Node was Fixed -------------------------------
		for (int i2=0; i2<fixQty; i2++) {
			if (fixEqu[i2] == i1) {
				f[i1] = 0;
				fixNode = 1;
			}
		}
		
		// --- Step A.3B: Creates Rates & Functions ------------------------------
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
 ** The parameters & variables needed for this function are declared in
 ** Function D (the Main Function) and are linked to this one. They include:
 ** t				double	The initial time
 ** t1			double	The final time (when the model should stop)
 ** maxTime		double	This should be the same value as 't1'
 ** maxNode		int		The maximum amount of nodes used.
 **							(note, maxNode is also declared seperately in Func A)
 ** y				double	The maxNode sized array which holds the final results
 **
 ** The Algorithm used for this function is:
 ** 1. We first declare the ODE system & the ODE driver which GSL requires.
 ** 2. We declare the parameters needed to save the output file.
 ** 3. A forLoop is declared to state how frequent the results should be saved.
 ** 3A. The stepsize is calculated
 ** 3B. The ODE solver is ran
 ** 3C. The results are displayed to the terminal window.
 ** 3D. The results are saved to the output file
 ** 4. The ODE solver is stopped & the output file is saved.
 ** ------------------------------------------------------------------------- */
void odePass (double t, double t1, double maxTime, double y[], int intRun,
				  int intMax) {
	
	// --- Step C.1: Declare ODE System & Driver --------------------------------
	gsl_odeiv2_system sys = {func, jac, maxNode};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new
		(&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);
	
	// --- Step C.2: Declare Output File Parameters -----------------------------
	char outBuffer[50];
	sprintf(outBuffer, "output_%d.dat", intRun);
	FILE *outFile = fopen(outBuffer, "w");
	outFile = fopen(outBuffer, "w");
	//fprintf(outFile, "Time\t Y0\t Y1\t Y2\t Y3\t Y4\t Y5\t Y6\t Y7\n");
	
	// --- Step C.3: States how often should Results be Saved -------------------
	for (int i=1; i<=maxTime; i++) {
      
		// --- Step C.3A: Calculate Stepsize -------------------------------------
		double ti = i * t1 / maxTime;
		
		// --- Step C.3B: ODE Solver ---------------------------------------------
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
		if (status != GSL_SUCCESS) {
			printf ("error, return value=%d\n", status);
			break;
		}
		
		// --- Step C.3C: Displays Results to Terminal ---------------------------
		//printf ("\nTime: %f\n",t);
		//printf ("%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
		//		  y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
		
		// --- Step C.3D: Saves Results to an Output File ------------------------
		fprintf (outFile,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",
					t, y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
	}
	
	
	// --- Step C.4: Ends ODE Solver & Saves Output File ------------------------
	gsl_odeiv2_driver_free (d);
	fclose(outFile);
	//fclose(sumFile);
	
}	// Function C Terminates Here



/* --- FUNCTION D: MAIN_PROGRAM ------------------------------------------------
 ** This program is the 'main' program of this script. This controls how the
 ** equations are modified before processing by the ODE solver.
 **
 ** These Func. C defined variables will be used: t, t1, maxTime, maxNode, & y
 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
 **
 ** The parameters & variables are consistent for all Func D Query are:
 ** intMethod	int		Saves the method the user selected in Step D.2.
 ** yInitial	double	A [maxNode] sized array with Y's initial values.
 ** y				double	A [maxNode] sized array with Y's initial values,
 **							after the user modifies its values.
 **
 ** The parameters & variables that will be modified for each Func D Query are:
 ** fixNode		int		States which nodes are fixed.
 ** fixVal		double	States a node's initial value.
 ** intRun		int		Tracks how often the ODE Solver (Func C) was ran.
 ** intMax		int		Predicts how many runs (intRun) will be required.
 **
 ** The Algorithm used for this function is:
 ** 1. We first declare the parameters & variables which Functions D & C need.
 ** 2. The user gets to select which query they'd like to run & how results
 **	 should be presented.
 **	 (Further documentation for thse queries are noted prior to each query)
 ** ------------------------------------------------------------------------- */
int main (int argc, char **argv) {

	// --- Step D.1: Declare Parameters & Variables -----------------------------
	int intMethod = 0;	// States the Method the user picked
	int intRun = 0;		// States the Run Number the program is on
	int intMax = 1;		// States the Max Amount of Runs Predicted
	
	double t = 0.0;									// Initial Time for Model
	double maxTime = 200, t1 = maxTime;			// Max Time for the Model
	
	double yInit[8] =							// Initial Conditions for Y
	{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	double y[maxNode];								// Creates Array for Modified Y
	
	for (int i=0; i<maxNode; i++)					// Sets Initial Y = Modified Y
		y[i] = yInit[i];
	
	// --- Step D.2: Selection of Query -----------------------------------------
	printf("Select a Method: "); scanf("%d",&intMethod);
	
	/* --- QUERY D.2.0: Manual ODE ----------------------------------------------
	 ** This query lets the user run the ODE Function & fix specified nodes.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 **
	 ** The following parameters & variables will be modified in for Query 0:
	 ** fixVal		User Set		States the value each node should be fixed with.
	 **
	 ** The following parameter & variables are unique to Query 0:
	 ** fixNode		User Set		States which nodes should be fixed.
	 **
	 ** The Algorithm used for Query 0 is:
	 ** 0A. The user gets to select how many nodes should be fixed
	 ** 0B. An ifStatement is declared if the user said a node should be fixed.
	 **	 In it, a forLoop is declared to have the user specify which nodes
	 **	 should be fixed and with what initial value.
	 ** 0C. The odeFunc (Func C) is ran to produce results.
	 ** ---------------------------------------------------------------------- */
	if (intMethod == 0) {
		// --- Step D.2.0A: Asks for Qty of Nodes to Fix -------------------------
		int fixNode = 0; double fixVal;
		printf("\nHow many nodes should be fixed? = "); scanf("%d",&fixQty);
		
		// --- Step D.2.0B: Asks for which Nodes to Fix & how so -----------------
		if (fixQty >= 1) {
			printf("\nThe range of nodes is 0 to %d",maxNode-1);
			for (int i=0; i<fixQty; i++) {
				printf("\nKnockout Node Number= "); scanf("%d",&fixNode);
				fixEqu[i] = fixNode;
				printf("Knockout with the (double) Value = "); scanf("%lf",&fixVal);
				y[fixNode] = fixVal;
				intRun++;
			}
		}
		
		// --- Step D.2.0C: Runs ODE Solver --------------------------------------
		odePass(t, t1, maxTime, y, intRun, intMax);
	}

	
	/* --- QUERY D.2.1: Automated SKO -------------------------------------------
	 ** This query lets the user do an automated Single Knock Out routine.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 **
	 ** The following parameters & variables will be modified in Query 1:
	 ** fixVal		Script Set		States the value a node will be fixed with.
	 **
	 ** ---------------------------------------------------------------------- */
	if (intMethod == 1) {
		fixQty = 1;								// SKO -- This fixes 1 node at a time
		intMax = maxNode * fixQty;
		for (int i1=0; i1<maxNode; i1++) {
			fixEqu[0] = i1; y[i1] = 0;		// Fixes the First Node Automatically
			odePass(t, t1, maxTime, y, intRun, intMax);	// Runs the ODE Solver
			printf("Node %d Done\n\n",fixEqu[0]);
			intRun++;
			y[i1] = yInit[i1];
		}
	}
	
	
	/* --- QUERY D.2.2: Automated DKO -------------------------------------------
	 ** This query lets the user do an automated Double Knock Out routine.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 **
	 ** The following parameters & variables will be modified in Query 2:
	 ** fixVal		Script Set		States the value a node will be fixed with.
	 **
	 ** ---------------------------------------------------------------------- */
	if (intMethod == 2) {
		fixQty = 2;								// DKO -- This fixes 2 nodes at a time
		intMax = pow(maxNode, fixQty);
		for (int i1=0; i1<maxNode; i1++) {
			fixEqu[0] = i1; y[i1] = 0;		// Fixes the First Node
			
			for (int i2=0; i2<maxNode; i2++) {
				fixEqu[1] = i2; y[i2] = 0;	// Fixes the Second Node
				odePass(t, t1, maxTime, y, intRun, intMax);	// Runs the ODE Solver
				printf("Node %d & %d Done\n\n",fixEqu[0],fixEqu[1]);
				intRun++;
				y[i2] = yInit[i2];			// Resets Node 2 to its original value
			}
			
			y[i1] = yInit[i1];				// Resets Node 1 to its Original Value
		}
	}
	

	/* --- QUERY D.2.4: Automated One-Node Sweep --------------------------------
	 ** This query lets the user do an automated sweep of one node.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 **
	 ** The following parameters & variables will be modified in Query 4:
	 ** fixVal		Script Set		States the value a node will be fixed with.
	 **
	 ** The following parameters & variables are unique for Query 4:
	 ** fixNode		User Set			States which nodes should be fixed.
	 ** maxVal		User Set			States the max value the fixed node will reach.
	 **
	 ** ---------------------------------------------------------------------- */
	if (intMethod == 4) {
		fixQty = 1;								// SKO -- This fixes 1 node at a time
		int fixNode = 0;
		double maxVal;
		
		printf("Fix Node # (Note: The 1st node is node 0) = ");
			scanf("%d",&fixNode);
		fixEqu[0] = fixNode;
		
		printf("Select the Max Value of the Fixed Node ");
			scanf("%lf",&maxVal);
		
		intMax = maxVal * (1/0.5);
		
		for (double i=0; i<=maxVal; i=i+0.1) {
			y[fixNode] = i;		// Fixes the First Node Automatically
			odePass(t, t1, maxTime, y, intRun, intMax);	// Runs the ODE Solver
			printf("Node %d at Value %f Done\n\n",fixEqu[0],i);
			intRun++;
			
		}
	}
	
	
	/* --- QUERY D.2.5: Automated Two-Node Sweep --------------------------------
	 ** This query lets the user do an automated sweep of two nodes.
	 **
	 ** These globally defined variables will be used: fixQty, maxNode, fixEqu
	 **
	 ** The following parameters & variables will be modified in Query 5:
	 ** fixVal		Script Set		States the value a node will be fixed with.
	 **
	 ** The following parameters & variables are unique for Query 5:
	 ** fixNode		User Set			States which nodes should be fixed.
	 ** tempVal		User Set			States the max value the fixed node will reach.
	 ** maxVal		Script Set		A [fixQty] sized array storing the max values
	 **									all nodes will reach.
	 **
	 ** ---------------------------------------------------------------------- */
	if (intMethod == 5) {
		fixQty = 2;								// DKO -- This fixes 2 nodes at a time
		int fixNode = 0;
		double tempVal, maxVal[2];
		
		for (int i=0; i<2; i++) {
			printf("Fix Node # (Note: The 1st node is node 0) = ");
			scanf("%d",&fixNode);
			fixEqu[i] = fixNode;
			printf("Select the Max Value of the Fixed Node = ");
			scanf("%lf",&tempVal);
			maxVal[i] = tempVal;
		}
		
		intMax = (maxVal[0] * maxVal[1]);
		
		for (double i1=0; i1<=maxVal[0]; i1++) {
			y[fixEqu[0]] = i1;		// Fixes the First Node Automatically
			printf("\nNode %d is at %f\n",fixEqu[0],i1);
			for (double i2=0; i2<=maxVal[1]; i2++) {
				y[fixEqu[1]] = i2;		// Fixes the Second Node Automatically
				printf("Node %d is at %f\n",fixEqu[1],i2);
				odePass(t, t1, maxTime, y, intRun, intMax);	// Runs the ODE Solver
				intRun++;
			}
		}
	}
	
}	// Function D Terminates Here