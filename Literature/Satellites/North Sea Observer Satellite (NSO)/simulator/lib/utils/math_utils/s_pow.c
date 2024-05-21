
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

#define S_FUNCTION_NAME s_pow
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
/* #define k (double)0 */

/*=====================*
 *Calculation functions*
 *=====================*/

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
  ssSetNumSFcnParams(S,0); /* Number of expected parameters */

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
    set number of inputs
  */
  if (!ssSetNumInputPorts(S, 2)) return;

  /*
   * Set direct feedthrough flag (1=yes, 0=no).
   * A port has direct feedthrough if the input is used in either
   * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
   * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
   */

  /* Port 0 a signal, u of 1 dimension */
  ssSetInputPortWidth(S,  0, 1);
  ssSetInputPortDataType(S, 0, SS_DOUBLE); 
  ssSetInputPortDirectFeedThrough(S, 0, 1);
  ssSetInputPortRequiredContiguous(S, 0, 1); /*direct input signal access*/
    
  /* Port 1 another signal, v of 1 dimension */
  ssSetInputPortWidth(S,  1, 1);
  ssSetInputPortDataType(S, 1, SS_DOUBLE);
  ssSetInputPortDirectFeedThrough(S, 1, 1);
  ssSetInputPortRequiredContiguous(S, 1, 1); /*direct input signal access*/

  /*
    set number of outputs
  */
  if (!ssSetNumOutputPorts(S, 1)) return;

  /* Port 0 the cross product of u and v */
  ssSetOutputPortWidth(S, 0, 1);
  ssSetOutputPortDataType(S, 0, SS_DOUBLE);

  ssSetNumSampleTimes(S, 1);
  ssSetNumRWork(S, 0);
  ssSetNumIWork(S, 0);
  ssSetNumPWork(S, 0);
  ssSetNumModes(S, 0);
  ssSetNumNonsampledZCs(S, 0);

  ssSetOptions(S, 0);

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
/*static void mdlStart(SimStruct *S)
{*/
  /*
    This function updates the parameter dependent konstants
  */
/*}*/

/* Function: mdlOutputs =======================================================
 * Abstract:
 *    In this function, you compute the outputs of your S-function
 *    block. Generally outputs are placed in the output vector, ssGetY(S).
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
	real_T       *u = (real_T *) ssGetInputPortSignal(S,0);
	real_T       *v = (real_T *) ssGetInputPortSignal(S,1);
	real_T       *y = ssGetOutputPortSignal(S,0);	
	
	/*
	  Calculate u^v
	*/

	y[0] = pow(u[0],v[0]);
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












