
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
#include "var_global.h"
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 2
#define u_1_width 1
#define u_2_width 1
#define u_3_width 1
#define u_4_width 1
#define u_5_width 1
#define u_6_width 1
#define y_width 5

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* === Variables globales persistentes === */
extern double P[5][5];
extern double Q[5][5];
extern double x_prev[5];
extern double x_pred[5];
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void ekf_predictive_Outputs_wrapper(const real_T *u,
			const real_T *R,
			const real_T *L,
			const real_T *Ke,
			const real_T *J,
			const real_T *Nr,
			const real_T *Ts,
			real_T *x_pred_out)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
/* ----------------------------------Predicción ------------------------------- */ 
    
    // Estados previos
    double Id_k = x_prev[0];
    double Iq_k = x_prev[1];
    double Wm_k = x_prev[2];
    double theta_m_k = x_prev[3];
    double Tx_k = x_prev[4];

    // Relaciones eléctricas
    double We_k = (*Nr) * Wm_k;
    double theta_e_k = (*Nr) * theta_m_k;

    double Vd = u[1];
    double Vq = u[2];


    // === Ecuaciones diferenciales ===
    double did_dt      = (Vd - (*R) * Id_k + We_k * (*L) * Iq_k) / (*L);
    double diq_dt      = (Vq - (*R) * Iq_k - We_k * (*L) * Id_k - (*Ke) * We_k) / (*L);
    double dwm_dt      = ((*Ke) * Iq_k - Tx_k) / (*J);
    double dtheta_m_dt = Wm_k;
    double dTx_dt      = 0.0;

    // === Predicción discreta (Euler hacia adelante) ===
    x_pred[0] = Id_k      + (*Ts) * did_dt;
    x_pred[1] = Iq_k      + (*Ts) * diq_dt;
    x_pred[2] = Wm_k      + (*Ts) * dwm_dt;
    x_pred[3] = theta_m_k + (*Ts) * dtheta_m_dt;
    x_pred[4] = Tx_k      + (*Ts) * dTx_dt;

    // === Jacobiano A ===
    double A[5][5] = {
        { -(*R)/(*L),           We_k,           0,              Iq_k,      0 },
        { -We_k,            -(*R)/(*L),     -(*Ke)/(*L),       -Id_k,      0 },
        { 0,                (*Ke)/(*J),         0,                0,     -1/(*J) },
        { 0,                    0,              0,                0,       0 },
        { 0,                    0,              0,                0,       0 }
    };

    // === Matriz identidad ===
    double I[5][5] = {
        {1,0,0,0,0},
        {0,1,0,0,0},
        {0,0,1,0,0},
        {0,0,0,1,0},
        {0,0,0,0,1}
    };

    // === Discretización: A_d = I + Ts*A ===
    double A_d[5][5];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            A_d[i][j] = I[i][j] + (*Ts) * A[i][j];
        }
    }

    // === Propagación de la covarianza: P = A_d * P * A_d^T + Q ===
    double temp[5][5];
    double P_pred[5][5];

    // temp = A_d * P
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            temp[i][j] = 0.0;
            for (int k = 0; k < 5; k++) {
                temp[i][j] += A_d[i][k] * P[k][j];
            }
        }
    }

    // P_pred = temp * A_d^T + Q
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            P_pred[i][j] = Q[i][j];
            for (int k = 0; k < 5; k++) {
                P_pred[i][j] += temp[i][k] * A_d[j][k];
            }
        }
    }

    // Actualiza la covarianza
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            P[i][j] = P_pred[i][j];
        }
    }

    // Guarda la predicción como nuevo estado previo
    for (int i = 0; i < 5; i++) {
        x_prev[i] = x_pred[i];
        x_pred_out[i] = x_pred[i];
    }
    
/* ----------------------------------Corrección ------------------------------- */
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


