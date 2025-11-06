
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
#define u_4_width 1
#define u_5_width 1
#define u_6_width 1
#define u_7_width 1
#define u_8_width 1
#define u_9_width 1
#define y_width 5
#define y_1_width 5

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define TWO_PI (2.0 * M_PI)

/* ---------- Utilidades ---------- */
static inline double wrap_2pi(double a) {
    a = fmod(a, TWO_PI);
    if (a < 0.0) a += TWO_PI;
    return a;
}

/* === Variables globales persistentes === */
static double x_prev[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   // [Id, Iq, Wm, theta_m, Tx]
static double x_pred[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
static double P[5][5] = {
    {1e-3, 0, 0, 0, 0},
    {0, 1e-3, 0, 0, 0},
    {0, 0, 1, 0, 0},
    {0, 0, 0, 1e-3, 0},
    {0, 0, 0, 0, 1e-3}
};
static double Q[5][5] = {
    {1e-6, 0, 0, 0, 0},
    {0, 1e-6, 0, 0, 0},
    {0, 0, 1e-3, 0, 0},
    {0, 0, 0, 1e-4, 0},
    {0, 0, 0, 0, 1e-9}
};
// --- Variables locales ---
   static double H[2][5] = {
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0}
    };

    // Ruido de medición (ajústalo según tus sensores)
   static double Rm[2][2] = {
        {1e-4, 0},
        {0, 1e-4}
    };
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Start function
 *
 */
void ekf_Start_wrapper(void)
{
/* %%%-SFUNWIZ_wrapper_Start_Changes_BEGIN --- EDIT HERE TO _END */
 
/* %%%-SFUNWIZ_wrapper_Start_Changes_END --- EDIT HERE TO _BEGIN */
}
/*
 * Output function
 *
 */
void ekf_Outputs_wrapper(const real_T *Va,
			const real_T *Vb,
			const real_T *R,
			const real_T *L,
			const real_T *Kt,
			const real_T *J,
			const real_T *Ts,
			const real_T *Ia_medido,
			const real_T *Ib_medido,
			const real_T *Nr,
			real_T *x_pred_out,
			real_T *x_corr_out)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
// === ESTADOS PREVIOS ===
double Ia_k = x_prev[0];
double Ib_k = x_prev[1];
double Wm_k = x_prev[2];
double th_m_k = x_prev[3];
double Tx_k = x_prev[4];

// === CONSTANTES Y VARIABLES AUXILIARES ===
double We_k  = (*Nr) * Wm_k;     // Velocidad eléctrica
double th_e_k = (*Nr) * th_m_k;  // Ángulo eléctrico

// === ECUACIONES DIFERENCIALES ===
// ¡Ojo! Corrijo: las funciones trigonométricas deben ir en minúsculas (sin “Th” mayúscula).
// También se usa Wm_k (no Wm), y todo debe multiplicarse con el operador * explícito.
double dia_dt = ((*Va) - (*R) * Ia_k + (*Kt) * (*Nr) * Wm_k * sin(th_e_k)) / (*L);
double dib_dt = ((*Vb) - (*R) * Ib_k - (*Kt) * (*Nr) * Wm_k * cos(th_e_k)) / (*L);
double dwm_dt = ((*Kt) * (Ib_k * cos(th_e_k) - Ia_k * sin(th_e_k)) - Tx_k) / (*J);
double dth_m_dt = Wm_k;
double dTx_dt = 0.0;

// === PREDICCIÓN DISCRETA (Euler hacia adelante) ===
x_pred[0] = Ia_k + (*Ts) * dia_dt;
x_pred[1] = Ib_k + (*Ts) * dib_dt;
x_pred[2] = Wm_k + (*Ts) * dwm_dt;
x_pred[3] = th_m_k + (*Ts) * dth_m_dt;
x_pred[4] = Tx_k + (*Ts) * dTx_dt;

// === JACOBIANO A ===
// Correcciones:
// - Todas las variables deben tener el operador * para los punteros.
// - Faltaban multiplicaciones y paréntesis.
// - El ángulo es th_e_k, no “Th”.
// - Kt, Nr, R, L, J deben desreferenciarse con (*).
// - Evita errores de signos y productos.

double A[5][5] = {
    { -(*R)/(*L),  0.0,   (*Kt)*(*Nr)*sin(th_e_k)/(*L),   (*Kt)*(*Nr)*Wm_k*cos(th_e_k)/(*L),   0.0 },
    {  0.0,       -(*R)/(*L),  -(*Kt)*(*Nr)*cos(th_e_k)/(*L),  (*Kt)*(*Nr)*Wm_k*sin(th_e_k)/(*L),  0.0 },
    { -(*Kt)*sin(th_e_k)/(*J),  -(*Kt)*cos(th_e_k)/(*J),   0.0,
      -(*Kt)*(*Nr)/(*J)*(Ia_k*cos(th_e_k)+Ib_k*sin(th_e_k)),  -1.0/(*J) },
    {  0.0,  0.0,  1.0,  0.0,  0.0 },
    {  0.0,  0.0,  0.0,  0.0,  0.0 }
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
            A_d[i][j] = I[i][j] + Ts[0] * A[i][j];
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
            S[i][j] = Rm[i][j];
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
    y[0] = Ia_medido[0] - Hx[0];
    y[1] = Ib_medido[0] - Hx[1];

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
            temp5[i][j] = I[i][j] - KH[i][j];
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
        x_corr_out[i] = x_corr[i];
    }

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            P[i][j] = P_new[i][j];
        }
    }
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


