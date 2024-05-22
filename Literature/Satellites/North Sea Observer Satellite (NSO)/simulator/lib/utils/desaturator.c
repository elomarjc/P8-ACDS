
/*
  Specification on the S-function name and level.
  this s-function is compatible with simulink3.
 */

#define S_FUNCTION_NAME desaturator
#define S_FUNCTION_LEVEL 2

#include "simstruc.h"

#include "math.h"

/*real_T umeas[3],umeasold[3],ytemp[3];*/
char onstate[3],offstate[3],outstate;
real_T ytemp[3],btemp[3];



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
    setup the right number of parameters: wc C per td
  */
  ssSetNumSFcnParams(S,5);

  if(ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)){
    /*
      sfunction is initialised with the wrong number of parameters
    */
    return;
  }

  /*
    set number of continnuert and discrete states.
  */
  ssSetNumContStates(S, 0); /*3 continuert states to implement the controller*/
  ssSetNumDiscStates(S, 0); /*zero discrete states*/


  /*
    specify the number of inputs to 1
  */
  if (!ssSetNumInputPorts(S, 2)) return; /*wrong number of inputs*/
  
  ssSetInputPortWidth(S, 0, 3); 
  ssSetInputPortRequiredContiguous(S, 0, true); /*direct input signal access*/
  ssSetInputPortWidth(S, 1, 3); 
  ssSetInputPortRequiredContiguous(S, 1, true); /*direct input signal access*/

  /*
   * Set direct feedthrough flag (1=yes, 0=no).
   * A port has direct feedthrough if the input is used in either
   * the mdlOutputs or mdlGetTimeOfNextVarHit functions.
   * See matlabroot/simulink/src/sfuntmpl_directfeed.txt.
   */
  ssSetInputPortDirectFeedThrough(S, 0, 1); /*Direct feed forward is not used*/
  ssSetInputPortDirectFeedThrough(S, 1, 1); /*Direct feed forward is not used*/


  /*
    check if the output is setup right:
  */
  if (!ssSetNumOutputPorts(S, 1)) return; /*no output set*/
  ssSetOutputPortWidth(S, 0, 3);
  ssSetOutputPortReusable(S,0,1); /*This states that the output is persistent, neede for the use of a merge of the output. NOTE reset this if merge is not used*/
  ssSetNumSampleTimes(S, 1); /* Set the sample time*/

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
    ssSetSampleTime(S, 0, 0.0);
    ssSetOffsetTime(S, 0, 0.0);

}

#define MDL_INITIALIZE_CONDITIONS
/* Function: mdlInitializeConditions ========================================
 * Abstract:
 *    Initialize the states
 */

static void mdlInitializeConditions(SimStruct *S){
	onstate[0] = 0;
	onstate[1] = 0;
	onstate[2] = 0;
	offstate[0] = 0;
	offstate[1] = 0;
	offstate[2] = 0;
	outstate = 0;
	ytemp[0] = 0;
	ytemp[1] = 0;
	ytemp[2] = 0;
	btemp[0] = 0;
	btemp[1] = 0;
	btemp[2] = 0;
}

/* Function: mdlOutputs =======================================================
 * Abstract:
 *      y = Cx + Du 
 */
static void mdlOutputs(SimStruct *S, int_T tid)
{
	real_T habs[3];

    const real_T *b    = (const real_T*)ssGetInputPortSignal(S,0);
    const real_T *h    = (const real_T*)ssGetInputPortSignal(S,1);
    real_T *y  = (real_T *)ssGetOutputPortSignal(S,0);

	real_T *on = (real_T *)mxGetPr(ssGetSFcnParam(S,1));
    real_T *off = (real_T *)mxGetPr(ssGetSFcnParam(S,2));
    real_T *tp = (real_T *)mxGetPr(ssGetSFcnParam(S,3));
    real_T *td = (real_T *)mxGetPr(ssGetSFcnParam(S,4));
    real_T *C = (real_T *)mxGetPr(ssGetSFcnParam(S,0));
    const real_T time =ssGetT(S);
    UNUSED_ARG(tid); /* not used in single tasking mode */

	habs[0] = h[0]>0 ? h[0] : -h[0];
	habs[1] = h[1]>0 ? h[1] : -h[1];
	habs[2] = h[2]>0 ? h[2] : -h[2];
	onstate[0] = habs[0]>*on;
	onstate[1] = habs[1]>*on;
	onstate[2] = habs[2]>*on;
	offstate[0] = habs[0]>*off;
	offstate[1] = habs[1]>*off;
	offstate[2] = habs[2]>*off;
	outstate = (outstate && (offstate[0]||offstate[1]||offstate[2]) ) || onstate[0]||onstate[1]||onstate[2];
	
	/*if(outstate) printf("outstate: %d\n",outstate);*/
	if(outstate) {
		if(fmod((double)time*10,*tp)<*td){
			btemp[0]=b[0];
			btemp[1]=b[1];
			btemp[2]=b[2];
			y[0]=0;
			y[1]=0;
			y[2]=0;
		}else{
			y[0]=-*C*(btemp[1]*h[2]-btemp[2]*h[1]);
			y[1]=-*C*(btemp[2]*h[0]-btemp[0]*h[2]);
			y[2]=-*C*(btemp[0]*h[1]-btemp[1]*h[0]);
		}
	}
	else {
		y[0]=0;
		y[1]=0;
		y[2]=0;
	}

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

















