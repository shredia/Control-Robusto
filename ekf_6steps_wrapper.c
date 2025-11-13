
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
#include <stdio.h>
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
static inline double wrap_2pi(double angle) {
    angle = fmod(angle + M_PI, TWO_PI);
    if (angle < 0.0) angle += TWO_PI;
    return angle - M_PI;
}


#define N 5

void matmul5x5(double A[N][N], double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
}
void transpose5x5(double A[N][N], double AT[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            AT[j][i] = A[i][j];
        }
    }
}

void matmul_f64(const double *A,
                const double *B,
                double *C,
                int m, int n, int p) //usamos las matrices C = A*B, A: mxn B:nxp C:mxp
{
    int i, j, k;
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            double sum = 0.0;
            for (k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * p + j];
            }
            C[i * p + j] = sum;
        }
    }
}


// Inversa de una matriz 2x2
// Entrada:  S[2][2]
// Salida:   Sinv[2][2]
// Si el determinante es 0, no se invierte (devuelve 0)
static inline int inv2x2(double S[2][2], double Sinv[2][2])
{
    double det = S[0][0]*S[1][1] - S[0][1]*S[1][0];
    if (det == 0.0) {
        return 0; // No se puede invertir
    }

    double invDet = 1.0 / det;

    Sinv[0][0] =  S[1][1] * invDet;
    Sinv[0][1] = -S[0][1] * invDet;
    Sinv[1][0] = -S[1][0] * invDet;
    Sinv[1][1] =  S[0][0] * invDet;

    return 1; // Éxito
}

/* === Variables globales persistentes === */
// [Id, Iq, Wm, theta_m, Tx]
static double x_hat_less[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   //estado predecido (con los k-1 estados)
static double x_dot_hat_less[5] = {0.0, 0.0, 0.0, 0.0, 0.0};//la derivada
static double x_hat_more[5] = {0.0, 0.0, 0.0, 0.0, 0.0};   //estado estimado  (con los k estados + entrada)
static double sum_less[5][5];
static double sum_more[5][5];
static double y_hat[2];
static double Lk[5][2];
static double y_niato[2];

static double Ck[2][5] = {
    {1.0,0.0,0.0,0.0,0.0 },
    {0.0,1.0,0.0,0.0,0.0}
};
static double CkT[5][2] = {
    {1.0,0.0},
    {0.0,1.0},
    {0.0,0.0},
    {0.0,0.0},
    {0.0,0.0}
};

static double Qk[5][5] = {
    {  1.0,  0.0,  0.0,  0.0,  0.0 },
    {  0.0,  1.0,  0.0,  0.0,  0.0 },
    {  0.0,  0.0,  1.0,  0.0,  0.0 },
    {  0.0,  0.0,  0.0,  1.0,  0.0 },
    {  0.0,  0.0,  0.0,  0.0,  1.0 }
};

//Covarianza del sensor, confiabilidad del sensor +- ,etc
static double sum_v[2][2] = { 
    { 1.0, 0.0},
    { 0.0, 1.0}
};
   
static double B = 0.01;
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output function
 *
 */
void ekf_6steps_Outputs_wrapper(const real_T *Va,
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
// === 1A: ESTADOS PREVIOS ===
double Ia_k =   x_hat_more[0];
double Ib_k =   x_hat_more[1];
double Wm_k =   x_hat_more[2];
double th_m_k = x_hat_more[3];
double Tx_k =   x_hat_more[4];

double th_e_k = wrap_2pi((*Nr) * th_m_k);
//calculamos las derivadas por EDO
//  x_hat_less (k) = f(x_hat_more (k-1),u(k-1),w(k-1))
x_dot_hat_less[0] = ((*Va) - (*R) * Ia_k + (*Kt)  * Wm_k * sin(th_e_k)) / (*L);
x_dot_hat_less[1] = ((*Vb) - (*R) * Ib_k - (*Kt)  * Wm_k * cos(th_e_k)) / (*L);
x_dot_hat_less[2] = ((*Kt)* (Ib_k * cos(th_e_k) - Ia_k * sin(th_e_k)) -B*Wm_k - Tx_k) / (*J);
x_dot_hat_less[3] = wrap_2pi(Wm_k);
x_dot_hat_less[4] = 0.0;
 //Integramos mediante euler 
x_hat_less[0] = x_hat_less[0] + (*Ts)*x_dot_hat_less[0]; 
x_hat_less[1] = x_hat_less[1] + (*Ts)*x_dot_hat_less[1];
x_hat_less[2] = x_hat_less[2] + (*Ts)*x_dot_hat_less[2];
x_hat_less[3] = x_hat_less[3] + (*Ts)*x_dot_hat_less[3];
x_hat_less[4] = x_hat_less[4] + (*Ts)*x_dot_hat_less[4];

     //Salida corregida
       x_pred_out[0] = x_hat_less[0];
       x_pred_out[1] = x_hat_less[1];
       x_pred_out[2] = x_hat_less[2];
       x_pred_out[3] = x_hat_less[3];
       x_pred_out[4] = x_hat_less[4];
// === 1B: Jacobiano ===
    //calculamos jacobiano en función de los estados 
    //x_niato_less (k) = x(k) - x_hat_less(k)
    //x(k) = serie de taylor
    //esto se transforma en x(k) = A_hat (k-1) x_niato_more (k-1) A_hat T (k-1) + B_hat (k-1) W_niato (k-1)
    //ahora, podemos encontrar la matriz de covarianza sum_less = A_hat*sum_more*A_hat T + B_hat*sum*B_hat T
    // Cómo no sabemos la relación entre B y Sum, simplemente definimos una matriz diagonal Qk, el cual representará la confianza en el modelo

    double A[5][5] = {
    { -(*R)/(*L),  0.0,   (*Kt)*sin(th_e_k)/(*L),   (*Kt)*(*Nr)*Wm_k*cos(th_e_k)/(*L),   0.0 },
    {  0.0,          -(*R)/(*L),  -(*Kt)*sin(th_e_k)/(*L),  -(*Kt)*(*Nr)*Wm_k*cos(th_e_k)/(*L),  0.0 },
    { -(*Kt)*sin(th_e_k)/(*J),  (*Kt)*cos(th_e_k)/(*J),   -B/(*J),
      -((*Kt)*(*Nr)/(*J))*(Ia_k*cos(th_e_k)+Ib_k*sin(th_e_k)),  -1.0/(*J) },
    {  0.0,  0.0,  1.0,  0.0,  0.0 },
    {  0.0,  0.0,  0.0,  0.0,  0.0 }
};

    // === Matriz identidad ===
    double I[5][5] = {
        {1.0,0.0,0.0,0.0,0.0},
        {0.0,1.0,0.0,0.0,0.0},
        {0.0,0.0,1.0,0.0,0.0},
        {0.0,0.0,0.0,1.0,0.0},
        {0.0,0.0,0.0,0.0,1.0}
    };

    // === Discretización: A_d = I + Ts*A ===
    double A_hat[5][5];
    double A_hat_T[5][5];
    double temp[5][5];
    
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            A_hat[i][j] = I[i][j] + Ts[0] * A[i][j];
             A_hat_T[j][i] = A_hat[i][j];
        }
    }
   

    // --- Multiplicaciones de matrices ---
// temp = A_hat * sum_more
matmul5x5(A_hat, sum_more, temp);

// sum_less = temp * A_hat_T
matmul5x5(temp, A_hat_T, sum_less);

// --- Sumar Qk --- cómo es solo diagonal 
for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
        sum_less[i][j] += Qk[i][j];
    }
}

// === 1C: predecir salida del sistema y (k) ===
    //aun que deberiamos tener una función y (k) = h(hk,uk,vk), ó de la forma: y (k) = C*x_hat_less + D*Uk, considerando que nuestra matriz D
    //es 0, y que C es solamente [1 0 0 0 0; 0 1 0 0 0] nos podemos ahorrar mucho cálculo al realizar esto:
//usamos las matrices C = A*B, A: mxn B:nxp C:mxp
    matmul_f64(&Ck[0][0], x_hat_less, y_hat,2,5,1);
    //Hay que tener en consideración que los siguientes pasos también se pueden ahorrar saltandonos H, pero puede provocar errores si no se hace bien


    // === 2A: Estimador de la matriz de ganancias de Kalman Lk ===
    //Lk es el punto más crítico, tiene demasiados cálculos. Es de la siguiente forma:
    // Lk = sum_less * CkT*inv[Ck*sum_less*CkT+sum_v]


    //calculamos en una variable temporal 
    //temp_1 =Ck*sum_less
    double temp_1[2][5];
    matmul_f64(&Ck[0][0],&sum_less[0][0],&temp_1[0][0],2,5,5);
    //temp_2 = temp_1*sum_less
    double temp_2[2][2];
    matmul_f64(&temp_1[0][0],&CkT[0][0],&temp_2[0][0],2,5,2);
    //ahora sumamos temp_2 + sum_v (covarianza de incertidumbre del sensor)
    double sum_y[2][2];
    for(int i = 0; i< 2; i++){
        for(int j=0;j<2;j++){
            sum_y[i][j] = temp_2[i][j] + sum_v[i][j];
        }
    }
    double temp_4[2][2]; //calculamos la inversa del parentesis
    if(!inv2x2(sum_y,temp_4)){
        printf("error, det igual a 0");
              }
    
    //multiplicamos sum_less con CkT
     double temp_5[5][2];
     //usamos las matrices C = A*B, A: mxn B:nxp C:mxp (m,n,p)
    matmul_f64(&sum_less[0][0], &CkT[0][0], &temp_5[0][0],5,5,2);
    //terminamos de calcular Lk = temp_5*temp_4
    matmul_f64(&temp_5[0][0],&temp_4[0][0],&Lk[0][0],5,2,2);
    
    // === 2b: predecir salida del sistema y (k) ===
    //Calculamos ahora x_hat_more = x_hat_less + Lk*(yk-y_hat), dónde yk es la medición real y y_hat la medición que debería de ser

    //Calculamos el error de la salida y_niato = y_meas - y_hat
    y_niato[0] = Ia_medido[0] - y_hat[0];
    y_niato[1] = Ib_medido[0] - y_hat[1];
    //Calculamos la multiplicacion de la matriz Lk*y_niato [5x2][2x1]
    double temp_6[5][1];
    matmul_f64(&Lk[0][0],y_niato,&temp_6[0][0],5,2,1);
    
    //sumamos finalmente x_hat_more = x_hat_less + Lk*(y_niato)
    x_hat_more[0] = x_hat_less[0] + temp_6[0][0];
    x_hat_more[1] = x_hat_less[1] + temp_6[1][0];
    x_hat_more[2] = x_hat_less[2] + temp_6[2][0];
    x_hat_more[3] = wrap_2pi(x_hat_less[3] + temp_6[3][0]);
    x_hat_more[4] = x_hat_less[4] + temp_6[4][0];   
     
    //Salida corregida
       x_corr_out[0] = x_hat_more[0];
       x_corr_out[1] = x_hat_more[1];
       x_corr_out[2] = x_hat_more[2];
       x_corr_out[3] = x_hat_more[3];
       x_corr_out[4] = x_hat_more[4];
    // === 2C: Arreglar covarianza ===
    double LkT[2][5]; //calculamos la traspuesta
LkT[0][0] = Lk[0][0];
LkT[0][1] = Lk[1][0];
LkT[0][2] = Lk[2][0];
LkT[0][3] = Lk[3][0];
LkT[0][4] = Lk[4][0];

LkT[1][0] = Lk[0][1];
LkT[1][1] = Lk[1][1];
LkT[1][2] = Lk[2][1];
LkT[1][3] = Lk[3][1];
LkT[1][4] = Lk[4][1];
    //Actualizamos la covarianza mediante sum_more = sum_less - Lk*sum_y*LkT;
    double temp_7[5][2];
    matmul_f64(&Lk[0][0],&sum_y[0][0],&temp_7[0][0],5,2,2); //[5x2][2x2] 
    double temp_8[5][5];
    matmul_f64(&temp_7[0][0],&LkT[0][0],&temp_8[0][0],5,2,5); //[5x2][2x5]
    //sumamos
    for(int i = 0; i< 5; i++){
        for(int j=0;j<5;j++){
            sum_more[i][j] = sum_less[i][j] - temp_8[i][j];
        }
    }
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


