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
#define GM (double)(3.99*pow(10,14))

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

double calcR2(double R, double theta) {
	double result;
	result = -((3*GM*pow(R0,2))/(2*pow(R,4)))*J2*(3*cos(theta) - 1);
	return result;
}

double calcR3(double R, double theta) {
	double result;
	result = -((20*GM*pow(R0,3))/(2*pow(R,5)))*J3*(pow(cos(theta),3) -
			(3/5)*cos(theta));
	return result;
}

double calcR4(double R, double theta) {
	double result;
	result = -((175*GM*pow(R0,4))/(8*pow(R,6)))*J4*(pow(cos(theta),4) -
			(6/7)*pow(cos(theta),2) + (3/35));
	return result;
}

double calcangle2(double R, double theta) {
	double result, temp;
	temp = -((3*GM*pow(R0,2))/(pow(R,3)))*J2;
	result = cos(theta)*sin(theta);
	return temp*result;
	
}

double calcangle3(double R, double theta) {
	double result, temp;
	temp = -((5*GM*pow(R0,3))/(2*pow(R,4)))*J3;
	result = (- 3*pow(cos(theta),2)*sin(theta) + (3/5) * sin(theta));
	return temp*result;
}

double calcangle4(double R, double theta) {
	double result, temp;
	temp = -((35*GM*pow(R0,4))/(8*pow(R,5)))*J4;
	result = (- 4*pow(cos(theta),3)*sin(theta) + (12/7)*cos(theta)*sin(theta));
	return temp*result;
}
*/

double dx(x, y, z){
	double result;
	result=-2*(x/(pow((1+(pow(x,2)+pow(y,2))/z),2)*z));
	return result;
}

double dy(x, y, z){
	double result;
	result=-2*(y/(pow((1+(pow(x,2)+pow(y,2))/z),2)*z));
	return result;
}

double dz(x, y, z){
	double result;
	result=(pow(x,2)+pow(y,2))/(pow((1+(pow(x,2)+pow(y,2))/z),2)*pow(z,2));
	return result;
}

void calcp(double x, double y, double z, double *p, double R, double R2) {
	p[1]=atan(sqrt(fabs((pow(x,2)+pow(y,2))/z)));
    p[2]=pow(cos(p[1]),4)-(6/7)*pow(cos(p[1]),2)+(3/35);
	p[3]=pow(cos(p[1]),3)-((3/5)*cos(p[1]));
	p[4]=((1/2)*(pow(R0,2)*J2*(3*pow(cos(p[1]),2)-1)/R2))
			+((5/2)*((pow(R0,3)*J3*p[3])/pow(R2,(3/2)))) +((35/8)*((pow(R0,4)*J4*p[2])/pow(R2,2)));	
}

double Ux(double x, double y, double z, double p1, double p2, double p3, double p4, double p5, double R, double R2, double Rz) {
	double result, tmp1, tmp2;
	tmp1=(-GM*p5*x)/(pow(R2,(3/2)));
	
	tmp2=(-pow(R0,2)*J2*(pow(cos(p2),2) -1)*x)/(pow(R2,2));
	tmp2+=-3*((pow(R0,2)*J2*cos(p2)*sin(p2)*p1*x)/(R2*Rz));
	tmp2+=-((15/2)*((pow(R0,3)*J3*p4*x)/pow(R2,(5/2))));
	tmp2+=(5/2)*(pow(R0,3)*J3*((-3*((pow(cos(p2),2)*sin(p2)*p1*x)/Rz))+(3/5)*((sin(p2)*p1*x)/Rz))/pow(R2,(3/2)));
	tmp2+=-(35/2)*((pow(R0,4)*J4*p3*x)/pow(R2,3));
	tmp2+=(35/8)*((pow(R0,4)*J4*(-4*((pow(cos(p2),3)*sin(p2)*p1*x)/Rz)+(12/7)*((cos(p2)*sin(p2)*p1*x)/Rz)))/pow(R2,2));
	
	result=tmp1+((GM*tmp2)/R);
	return result;
}

double Uy(double x, double y, double z, double p1, double p2, double p3, double p4, double p5, double R, double R2, double Rz) {
	double result, tmp1, tmp2;
	tmp1=-((GM*p5*y)/(pow(R2,(3/2))));
	
	tmp2=(-pow(R0,2)*J2*(pow(cos(p2),2) -1)*y)/(pow(R2,2));
	tmp2+=-3*((pow(R0,2)*J2*cos(p2)*sin(p2)*p1*y)/(R2*Rz));
	tmp2+=-((15/2)*((pow(R0,3)*J3*p4*y)/pow(R2,(5/2))));
	tmp2+=(5/2)*(pow(R0,3)*J3*((-3*((pow(cos(p2),2)*sin(p2)*p1*y)/Rz))+(3/5)*((sin(p2)*p1*y)/Rz))/pow(R2,(3/2)));
	tmp2+=-(35/2)*((pow(R0,4)*J4*p3*y)/pow(R2,3));
	tmp2+=(35/8)*((pow(R0,4)*J4*(-4*((pow(cos(p2),3)*sin(p2)*p1*y)/Rz)+(12/7)*((cos(p2)*sin(p2)*p1*y)/Rz)))/pow(R2,2));
	
	result=tmp1+((GM*tmp2)/R);
	return result;
	
}

double Uz(double x, double y, double z, double p1, double p2, double p3, double p4, double p5, double R, double R2, double Rz) {
	double result, tmp1, tmp2, xy;
	tmp1=-((GM*p5*z)/(pow(R2,(3/2))));
	xy=pow(x,2)+pow(y,2);
	Rz=Rz*z;
	
	tmp2=(-pow(R0,2)*J2*(pow(cos(p2),2) -1)*z)/(pow(R2,2));
	tmp2+=(3/2)*((pow(R0,2)*J2*cos(p2)*sin(p2)*p1*xy)/(R2*Rz));
	tmp2+=-((15/2)*((pow(R0,3)*J3*p4*z)/pow(R2,(5/2))));
	tmp2+=(5/2)*(pow(R0,3)*J3*(((3/2)*((pow(cos(p2),2)*sin(p2)*p1*xy)/Rz))-(3/10)*((sin(p2)*p1*xy)/Rz)))/pow(R2,(3/2));
	tmp2+=-(35/2)*((pow(R0,4)*J4*p3*z)/pow(R2,3));
	tmp2+=(35/8)*((pow(R0,4)*J4*(2*((pow(cos(p2),3)*sin(p2)*p1*xy)/Rz)-(6/7)*((cos(p2)*sin(p2)*p1*xy)/Rz)))/pow(R2,2));
	
	result=tmp1+((GM*tmp2)/R);
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
	double x, y, z;
	double p[5];
	double R, R2, Rz;

	x=u[0];
	y=u[1];
	z=u[2];

	if(abs(z)<0.01){
	z=0.01;
	}	
	
	R=sqrt(fabs(pow(x,2)+pow(y,2)+pow(z,2)));
	R2=pow(x,2)+pow(y,2)+pow(z,2);
	Rz=sqrt(fabs((pow(x,2)+pow(y,2)/z)))*z;
			
	calcp(x, y, z, (double *)&p, R, R2);
	
	p[0]=dx(x,y,z);	
	result[0]=Ux(x, y, z, p[0], p[1], p[2], p[3], p[4], R, R2, Rz);
	p[0]=dy(x,y,z);	
	result[1]=Uy(x, y, z, p[0], p[1], p[2], p[3], p[4], R, R2, Rz);
	p[0]=dz(x,y,z);	
	result[2]=Uz(x, y, z, p[0], p[1], p[2], p[3], p[4], R, R2, Rz);

	
	
	
	/*sumr = calcR2(u[0],u[1]);
	sumr += calcR3(u[0],u[1]);
	sumr += calcR4(u[0],u[1]);
	
	sumangle = calcangle2(u[0],u[1]);
	sumangle += calcangle3(u[0],u[1]);
	sumangle += calcangle4(u[0],u[1]);
	
	result[0] = sumr;
	result[1] = 0;
	result[2] = sumangle;
	
	sum = calc2(u[0],u[1]);
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












