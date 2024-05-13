
/**
   This file implements Zonal harmonics in earths gravitational fields

   developed by:
      Kresten Kjeldgaard(kkje01@control.auc.dk)

   For the AAUSAT-II ADCS project.
 **/

/*
  Specification on the S-function name and level.
  this s-function is compatible with simulink3.
 */

#define S_FUNCTION_NAME zonal
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*====================*
 * IGRF model methods *
 *====================*/

/*
  constants
 */
#define J2 (double)0.00108263
#define J3 (double)-0.00000254
#define J4 (double)-0.00000161

#define R0 6370000
#define GM (double)3.99*pow(10,14)

/*=====================*
 *Calculation functions*
 *=====================*/
/*
double calc2(double R, double theta) {
	double result, temp;
	temp = (double)R0/R;
	result = (pow(temp,2))*J2;
	temp = .5*(3*(pow(cos(theta),2)) - 1);
    return result*temp;
}

double calc3(double R, double theta) {
	double result, temp;
	temp = (double)R0/R;
	result = (pow(temp,3))*J3;
	temp = (5/2)*(pow(cos(theta),3) - ((3/5)*cos(theta)));
    return result*temp;
}

double calc4(double R, double theta) {
	double result, temp;
	temp = (double)R0/R;
	result = (pow(temp,4))*J4;
	temp = (35/8)*(pow(cos(theta),4) - ((6/7)*pow(cos(theta),2)) + (3/35));
    return result*temp;
}
*/
double calcR12(double R, double theta) {
	double result;
	result = (1/2)*(((pow(R0,2))*J2*(3*pow(cos(theta),2) - 1))/pow(R,2));
	return result;
}

double calcR13(double R, double theta) {
	double result;
	result = (5/2)*(((pow(R0,3))*J3*(pow(cos(theta),3) - (3/5)*cos(theta)))/pow(R,3));
	return result;
}

double calcR14(double R, double theta) {
	double result;
	result = (35/8)*(((pow(R0,4))*J4*(pow(cos(theta),4) - (6/7)*pow(cos(theta),2) + (3/35)))/pow(R,4));
	return result;
}

double calcR22(double R, double theta) {
	double result;
	result = -(((pow(R0,2))*J2*(3*pow(cos(theta),2) - 1))/pow(R,3));
	return result;
}

double calcR23(double R, double theta) {
	double result;
	result = -(15/2)*(((pow(R0,3))*J3*(pow(cos(theta),3) - (3/5)*cos(theta)))/pow(R,4));
	return result;
}

double calcR24(double R, double theta) {
	double result;
	result = -(35/2)*(((pow(R0,4))*J4*(pow(cos(theta),4) - (6/7)*pow(cos(theta),2) + (3/35)))/pow(R,5));
	return result;
}

double calcangle2(double R, double theta) {
	double result;
	result = -3*(((pow(R0,2))*J2*(cos(theta)*sin(theta)))/pow(R,2));
	return result;
	
}

double calcangle3(double R, double theta) {
	double result;
	result = (5/2)*((pow(R0,3))*J3*(-3*pow(cos(theta),2)*sin(theta) + (3/5)*sin(theta))/pow(R,3));
	return result;
}

double calcangle4(double R, double theta) {
	double result;
	result = (35/8)*(((pow(R0,4))*J4*(-4*pow(cos(theta),3)*sin(theta) + (12/7)*cos(theta)*sin(theta)))/pow(R,4));
	return result;
}


/*====================*
 * S-function methods *
 *====================*/

/* Function: mdlInitializeSizes ===============================================
 * Abstract:
 *    The sizes information is used by Simulink to determine the S-function
 *    block's characteristics (number of inputs, outputs, states, etc.).
 */
static void mdlInitializeSizes(SimStruct *S)
{
  /*
    setup the right number of parameters: distance, coelevation
  */
  ssSetNumSFcnParams(S,0);

  if(ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)){
    /*
      sfunction is initialised with the wrong number of parameters
    */
    return;
  }

  /*
    set number of continnuert and discrete states.
  */
  ssSetNumContStates(S, 0); /*zero continuert states*/
  ssSetNumDiscStates(S, 0); /*zero discrete states*/


  /*
    specify the number of inputs to 1
    (|R_sc(I)| and coelecation)
  */
  if (!ssSetNumInputPorts(S, 1)) return; /*wrong number of inputs*/
  ssSetInputPortWidth(S, 0, 3);
  ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/

  /*
   * Set direct feedthrough flag (1=yes, 0=no).
   * A port has direct feedthrough if the input is used in either
   * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
   * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
   */
  ssSetInputPortDirectFeedThrough(S, 0, 1);

  /*
    check if the output is setup right:
    gravitional potential of earth zonal harmonics.
  */
  if (!ssSetNumOutputPorts(S, 1)) return; /*no output set*/
  ssSetOutputPortWidth(S, 0, 3);

  ssSetNumSampleTimes(S, 1);

}

/* Function: mdlInitializeSampleTimes =========================================
 * Abstract:
 *    This function is used to specify the sample time(s) for your
 *    S-function. You must register the same number of sample times as
 *    specified in ssSetNumSampleTimes.
 */
static void mdlInitializeSampleTimes(SimStruct *S)
{
  /*
    Set the sample time to continuous
  */
  ssSetSampleTime(S, 0, CONTINUOUS_SAMPLE_TIME);
  ssSetOffsetTime(S, 0, 0.0);
}


/* Function: mdlStart =======================================================
 * Abstract:
 *    This function is called once at start of model execution. If you
 *    have states that should be initialized once, this is the place
 *    to do it.
 */
static void mdlStart(SimStruct *S)
{
  /*
    This function updates the parameter dependent konstants
  */
}

/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block. Generally outputs are placed in the output vector, ssGetY(S).
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
	real_T       *u = (real_T *) ssGetInputPortSignal(S,0);
    real_T       *result = ssGetOutputPortSignal(S,0);	
	double sumr, sumr1, sumr2, sumtheta, sumangle, resultr, resultangle;
	
	sumr1 = calcR12(u[0],u[2]);
	sumr1 += calcR13(u[0],u[2]);
	sumr1 += calcR14(u[0],u[2]);
	sumr1 = -(GM*sumr1)/(pow(u[0],2));

	sumr2 = calcR22(u[0],u[2]);
	sumr2 += calcR23(u[0],u[2]);
	sumr2 += calcR24(u[0],u[2]);
	sumr2 = (GM*sumr2)/u[0];

	sumr = sumr1+sumr2;
	
	sumtheta = calcangle2(u[0],u[2]);
	sumtheta += calcangle3(u[0],u[2]);
	sumtheta += calcangle4(u[0],u[2]);
	
	sumangle = (GM*sumtheta)/(pow(u[0],2)*sin(u[1]));
	
	result[0] = sumr;
	result[1] = 0;
	result[2] = sumangle;
	
	/*sum = calc2(u[0],u[1]);
	sum += calc3(u[0],u[1]);
	sum += calc4(u[0],u[1]);

	y[0] = ((double)GM/u[0])*sum;*/	
	
			
	/*calc sum(7,u,B);
    *norm=l2norm(B,3);
    y[0]=B[0];
    y[1]=B[1];
    y[2]=B[2];*/
	
}


/* Function: mdlTerminate =====================================================
 * Abstract:
 *    In this function, you should perform any actions that are necessary
 *    at the termination of a simulation.  For example, if memory was
 *    allocated in mdlStart, this is the place to free it.
 */
static void mdlTerminate(SimStruct *S)
{
  /*
    Here the memory, used to store temporary data for the model is deallocated.
  */


}


/*=============================*
 * Required S-function trailer *
 *=============================*/

#ifdef  MATLAB_MEX_FILE    /* Is this file being compiled as a MEX-file? */
#include "simulink.c"      /* MEX-file interface mechanism */
#else
#include "cg_sfun.h"       /* Code generation registration function */
#endif












