
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
#define y_width 1

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
extern double P[5][5];
extern double Q[5][5];
extern double x_prev[5];
extern double x_pred[5];
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void ekf_correccion_Outputs_wrapper(const real_T *id_meas,
			const real_T *iq_meas,
			real_T *y0)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
// --- Variables locales ---
    double H[2][5] = {
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0}
    };

    // Ruido de medición (ajústalo según tus sensores)
    double R[2][2] = {
        {1e-4, 0},
        {0, 1e-4}
    };

    double I5[5][5] = {
        {1,0,0,0,0},
        {0,1,0,0,0},
        {0,0,1,0,0},
        {0,0,0,1,0},
        {0,0,0,0,1}
    };

    double Ht[5][2];          // Transpuesta de H
    double PHt[5][2];         // P * H^T
    double S[2][2];           // H * P * H^T + R
    double S_inv[2][2];       // Inversa de S
    double K[5][2];           // Ganancia de Kalman
    double y[2];              // Innovación (error de medición)
    double Hx[2];             // H * x_pred
    double KH[5][5];          // K * H
    double temp5[5][5];       // Matriz temporal

    // --- 1. Calcular H^T ---
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 2; j++) {
            Ht[i][j] = H[j][i];
        }
    }

    // --- 2. PHt = P * H^T ---
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 2; j++) {
            PHt[i][j] = 0.0;
            for (int k = 0; k < 5; k++) {
                PHt[i][j] += P[i][k] * Ht[k][j];
            }
        }
    }

    // --- 3. S = H * P * H^T + R ---
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            S[i][j] = R[i][j];
            for (int k = 0; k < 5; k++) {
                S[i][j] += H[i][k] * PHt[k][j];
            }
        }
    }

    // --- 4. Inversa de S (2x2) ---
    double detS = S[0][0]*S[1][1] - S[0][1]*S[1][0];
    S_inv[0][0] =  S[1][1] / detS;
    S_inv[0][1] = -S[0][1] / detS;
    S_inv[1][0] = -S[1][0] / detS;
    S_inv[1][1] =  S[0][0] / detS;

    // --- 5. K = P * H^T * S_inv ---
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 2; j++) {
            K[i][j] = 0.0;
            for (int k = 0; k < 2; k++) {
                K[i][j] += PHt[i][k] * S_inv[k][j];
            }
        }
    }

    // --- 6. y = z - H * x_pred ---
    Hx[0] = 0.0;
    Hx[1] = 0.0;
    for (int j = 0; j < 5; j++) {
        Hx[0] += H[0][j] * x_pred[j];
        Hx[1] += H[1][j] * x_pred[j];
    }
    y[0] = (*id_meas) - Hx[0];
    y[1] = (*iq_meas) - Hx[1];

    // --- 7. x_corr = x_pred + K * y ---
    double x_corr[5];
    for (int i = 0; i < 5; i++) {
        x_corr[i] = x_pred[i];
        for (int j = 0; j < 2; j++) {
            x_corr[i] += K[i][j] * y[j];
        }
    }

    // --- 8. P = (I - K * H) * P ---
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            KH[i][j] = 0.0;
            for (int k = 0; k < 2; k++) {
                KH[i][j] += K[i][k] * H[k][j];
            }
        }
    }

    // temp5 = I - K*H
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            temp5[i][j] = I5[i][j] - KH[i][j];
        }
    }

    // P = (I - K*H) * P
    double P_new[5][5];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            P_new[i][j] = 0.0;
            for (int k = 0; k < 5; k++) {
                P_new[i][j] += temp5[i][k] * P[k][j];
            }
        }
    }

    // --- Actualiza estados globales ---
    for (int i = 0; i < 5; i++) {
        x_prev[i] = x_corr[i];
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            P[i][j] = P_new[i][j];
        }
    }
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


