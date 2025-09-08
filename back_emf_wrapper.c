
/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#define SIMPLIFIED_RTWTYPES_COMPATIBILITY
#include "rtwtypes.h"
#undef SIMPLIFIED_RTWTYPES_COMPATIBILITY
#endif



/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <math.h>
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define u_1_width 1
#define u_2_width 1
#define u_3_width 1
#define y_width 1
#define y_1_width 1
#define y_2_width 1

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


static double R = 3.54;        // Ohm
static double L = 6.0/1000.0;  // H
static double Ts = 1e-3;       // s
static double Ialpha_prev = 0.0, Ibetha_prev = 0.0;
static double theta_prev = 0.0, we_prev = 0.0;
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void back_emf_Outputs_wrapper(const real_T *Valpha,
			const real_T *Vbetha,
			const real_T *Ialpha,
			const real_T *Ibetha,
			real_T *theta_e,
			real_T *w_e,
			real_T *dot_w_e)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
/* Derivadas de corriente */
double dIalpha_dt = (Ialpha[0] - Ialpha_prev) / Ts;
double dIbetha_dt = (Ibetha[0] - Ibetha_prev) / Ts;

/* Back-EMF (convención: Va = R Ia + L dIa/dt + Ea) */
double Ealpha = Valpha[0] - (R * Ialpha[0] + L * dIalpha_dt);
double Ebetha = Vbetha[0] - (R * Ibetha[0] + L * dIbetha_dt);

/* Ángulo eléctrico */
double theta = atan2(Ebetha, Ealpha);
*theta_e = theta;

/* Unwrap de ángulo para derivar correctamente */
double dtheta = theta - theta_prev;
if (dtheta >  M_PI) dtheta -= 2.0*M_PI;
if (dtheta < -M_PI) dtheta += 2.0*M_PI;

/* Velocidad eléctrica (rad/s) y su derivada */
double we = dtheta / Ts;
double dwe = (we - we_prev) / Ts;
*w_e      = we;
*dot_w_e  = dwe;

/* Actualizar memorias */
Ialpha_prev = Ialpha[0];
Ibetha_prev = Ibetha[0];
theta_prev  = theta;
we_prev     = we;


/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


