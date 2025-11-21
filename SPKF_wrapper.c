
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
#define y_width 4
#define y_1_width 4

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
// Necesitamos Na variables de estado para el sigma-point, que equivaldría a:
// NA = N + N_wk + N_Vk, dónde N son las variables de estado del sistema, N_vk la cantidad de variables ruidosas y N_vk la cantidad de variables que leemos
//en nuestro caso, Na = 4 + 4 + 2

#define NX 4
#define NW 4
#define NV 2
#define NA (NX+NW+NV)
#define N (2*NA+1)



/* === Variables globales persistentes === */
//X_hat_a es el vector de estados extendidos, compuesto por x_hat_more, w y v, siendo w y v las variables de estado extendidos.
//Aun que, en teoríam es un vector, la verdad es que es un conjunto de vectores, de tamaño 2n+1, dónde "barre" desde -n hasta n de la forma x +- y*sqrt(sum_a)
//sum_a es el conjunto de sigma-points generados del estado anterior corregido    

//Por lo tanto, X_a es un vector de tamaño Na x N (2*Na+1)
//para componer esta matriz, hay que componer los estados
static double x_hat[NX];      // estado estimado actual
static double x_hat_pred[NX]; // estado predicho
static double X_a[NA][N];        // sigma points actuales
static double X_a_pred[NA][N];   // sigma points propagados
//matriz unica de X_a
static double x_a[NA];


//Para componer sum_a, hay que tener las covarianzas
//Esta matriz, contra más pequeño, confía más en el modelo
static double sum_w[NX][NX] = {
    {1e-4,0.0,0.0,0.0},
    {0.0,1e-4,0.0,0.0},
    {0.0,0.0,1e-8,0.0},
    {0.0,0.0,0.0,1e-6}
};

//Matriz de ruido de medición
//contra más grande, menos cree en el sensor
static double sum_v[NV][NV] = {
    {1e-6,0.0},
    {0.0,1e-6}
};


//para el factor de Cholesky, necesitamos una matriz P y una matriz Sa
//P es la matriz que relaciona el error de los estados: P = E[(x-x_hat)(x-x_hat)T]
static double P[NX][NX]; //matriz de covarianza
//Para poder construir los sigma-point, es necesario descomponer la matriz P en los factores de cholesky.
static double Sa[NA][NA]; //Matriz de cholesky(Sa*SaT = P)
static double Pa[NA][NA]; //Matriz de covarianza aumentada    
static double Pa_copy[NA][NA];//copia de la matriz de covarianza
static double P_pred[NX][NX];
static double dx[NX];
static double Y_a_pred[NV][N]; //matriz de mediciones
static double y_hat_pred[NV];
static double S_cov[NV][NV]; //covarianzas de las mediciones
static double dy[NV];
static double P_xy[NX][NV];//covarianza cruzada 
static double innovation[NV];
static double S_inv[2][2];
static double K_gain[NX][NV];
static double KS[NX][NV];


int i = 0;
int j = 0;
int k = 0;

//Variables para derivadas, integrales y estados del sistema
static double dIa = 0;
static double dIb = 0;
static double dWm = 0;
static double dTh_m = 0;
static double S = 0;
static double C = 0;
static double Ia = 0;
static double Ib = 0;
static double Wm = 0;
static double Th_m = 0;
static double Th_e = 0;

static double B = 0.005;
static double gamma = 0;
//ruido
static double w0 = 0;
static double w1 = 0;
static double w2 = 0;
static double w3 = 0;

static double Wmedia[N];//Media
static double Wcov[N];//covarianza
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Start function
 *
 */
void SPKF_Start_wrapper(void)
{
/* %%%-SFUNWIZ_wrapper_Start_Changes_BEGIN --- EDIT HERE TO _END */
double alpha = 1e-4;
    double beta  = 2.0;
    double kappa = 1.0;

    double lambda = alpha*alpha * (NA + kappa) - NA;
    gamma = sqrt(NA + lambda);   // <-- este gamma usa el static double de arriba

    // Pesos de la media y la covarianza
    Wmedia[0] = lambda / (NA + lambda);
    Wcov[0]   = Wmedia[0] + (1.0 - alpha*alpha + beta);

    for (int kk = 1; kk < N; kk++) {
        Wmedia[kk] = 1.0 / (2.0 * (NA + lambda));
        Wcov[kk]   = Wmedia[kk];
    }

    // Estado inicial y P inicial
    for (int i = 0; i < NX; i++) {
        x_hat[i] = 0.0;    // o tu condición inicial de corrientes, velocidad, ángulo
        P[i][j] = 0.0;
    }
    P[0][0] = 1e-3;
    P[1][1] = 1e-3;
    P[2][2] = 1e-4;
    P[3][3] = 1e-6;

   /***** INICIALIZACIÓN COMPLETA DE TODAS LAS MATRICES IMPORTANTES ******/

// Matrices NA x NA
for (i = 0; i < NA; i++) {
    for (j = 0; j < NA; j++) {
        Sa[i][j]      = 0.0;
        Pa[i][j]      = 0.0;
        Pa_copy[i][j] = 0.0;
    }
}

// Matrices NX x NX
for (i = 0; i < NX; i++) {
    dx[i] = 0.0;
    for (j = 0; j < NX; j++) {
        P_pred[i][j] = 0.0;
    }
}

// Matrices NV x N
for (i = 0; i < NV; i++) {
    y_hat_pred[i] = 0.0;
    dy[i]         = 0.0;
    innovation[i] = 0.0;

    for (j = 0; j < N; j++) {
        Y_a_pred[i][j] = 0.0;
    
    }

    // Matrices NV x NV
    for (j = 0; j < NV; j++) {
        S_cov[i][j] = 0.0;
    }
}

// Matriz S_inv (2 x 2)
for (i = 0; i < 2; i++)
    for (j = 0; j < 2; j++)
        S_inv[i][j] = 0.0;

// Matrices NX x NV
for (i = 0; i < NX; i++) {
    for (j = 0; j < NV; j++) {
        P_xy[i][j]   = 0.0;
        K_gain[i][j] = 0.0;
        KS[i][j]     = 0.0;
    }
}

// Matrices NA x N (sigma points y predicciones)
for (i = 0; i < NA; i++) {
    for (j = 0; j < N; j++) {
        X_a[i][j]      = 0.0;
        X_a_pred[i][j] = 0.0;
    }
}

// Vector aumentado x_a
for (i = 0; i < NA; i++)
    x_a[i] = 0.0;
/* %%%-SFUNWIZ_wrapper_Start_Changes_END --- EDIT HERE TO _BEGIN */
}
/*
 * Output function
 *
 */
void SPKF_Outputs_wrapper(const real_T *Va,
			const real_T *Vb,
			const real_T *R,
			const real_T *L,
			const real_T *Kt,
			const real_T *J,
			const real_T *Ts,
			const real_T *Ia_medido,
			const real_T *Ib_medido,
			const real_T *Nr,
			real_T *x,
			real_T *x_aux)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
//1a-i: Calculamos X_hat_a en función de X_hat_more (estados)
    for(i = 0; i < NX; i ++){
        x_a[i] = x_hat[i];  
    }
    //calculamos en función de w (es 0) (ruido proceso)
    for(i = 0; i < NW; i ++){
        x_a[NX+i] = 0.0;  
    }
    //calculamos en función de v (es 0) (ruido medicion)
    for(i = 0; i < NV; i ++){
        x_a[NX+NW+i] = 0.0;  
    }
    // === 1A-ii: Cholesky de Pa  Sa ===


    // Construir Pa a partir de P, sum_w y sum_v
// Pa = [ P      0      0
//        0    sum_w   0
//        0     0    sum_v ]
// poner toda Pa en 0
for (i = 0; i < NA; i++)
  for (j = 0; j < NA; j++)
    Pa[i][j] = 0.0;
// Bloque P
for (i = 0; i < NX; i++)
    for (j = 0; j < NX; j++)
        Pa[i][j] = P[i][j];

// Bloque sum_w
for (i = 0; i < NW; i++)
    for (j = 0; j < NW; j++)
        Pa[NX + i][NX + j] = sum_w[i][j];

// Bloque sum_v
for (i = 0; i < NV; i++)
    for (j = 0; j < NV; j++)
        Pa[NX + NW + i][NX + NW + j] = sum_v[i][j];



    //construimos copia de Pa
for (int i = 0; i < NA; i++)
    for (int j = 0; j < NA; j++)
        Pa_copy[i][j] = Pa[i][j];
    



    //limpiamos Sa total
    for (int i = 0; i < NA; i++)
        for (int j = 0; j < NA; j++)
            Sa[i][j] = 0.0;
        
  


// Cholesky lower-triangular en Sa
for (int i = 0; i < NA; i++) {
    for (int j = 0; j <= i; j++) {

        double sum = Pa_copy[i][j];

        for (int k = 0; k < j; k++)
            sum -= Sa[i][k] * Sa[j][k];

        if (i == j) {
            if (sum <= 0.0) sum = 1e-12;
            Sa[i][j] = sqrt(sum);
        } else {
            if (fabs(Sa[j][j]) < 1e-15) {
    Sa[j][j] = 1e-6;  // parche de estabilidad
}
            Sa[i][j] = sum / Sa[j][j];
        }
    }

    //limpieza diagonal superior de la matriz
    for (int j = i+1; j < NA; j++)
        Sa[i][j] = 0.0; 
}

     // === 1A-iii: Construir los sigma-points X_a(:,k) a partir de x_a y Sa ===

// sigma point central
for (i = 0; i < NA; i++) {
    X_a[i][0] = x_a[i];
}

// sigma points positivos y negativos
for (j = 0; j < NA; j++) {
    for (i = 0; i < NA; i++) {
        X_a[i][1 + j]      = x_a[i] + gamma * Sa[i][j];
        X_a[i][1 + NA + j] = x_a[i] - gamma * Sa[i][j];
    }
}

    
    //Ahora aplicamos el ciclo principal, dónde produciremos 2n+1 iteraciones y propagamos el sigma point por el modelo
    for(k = 0; k<N;k++){
    //conseguimos las variables de estado
    Ia  = X_a[0][k];
    Ib  = X_a[1][k];
    Wm  = X_a[2][k];
    Th_m = X_a[3][k];
    Th_e = Th_m*(*Nr);

    w0 = X_a[NX + 0][k];
    w1 = X_a[NX + 1][k];
    w2 = X_a[NX + 2][k];
    w3 = X_a[NX + 3][k];

    //calculamos las derivadas sumandole el ruido
    S = sin((Th_e));
    C = cos((Th_e));
    dIa = (((*Va) - (*R)*(Ia) + (Wm)*(*Kt)*S)/(*L))+w0;
    dIb = (((*Vb) - (*R)*(Ib) - (Wm)*(*Kt)*C)/(*L))+w1;
    dWm = (((*Kt)*(Ib*C-Ia*S)-B*(Wm))/(*J))+w2;
    dTh_m = Wm;
    //integramos mediante el método de euler
     X_a_pred[0][k] = Ia  + (*Ts)*dIa;
     X_a_pred[1][k] = Ib  + (*Ts)*dIb;
     X_a_pred[2][k] = Wm  + (*Ts)*dWm;
     X_a_pred[3][k] = Th_m  + (*Ts)*dTh_m;
     X_a_pred[3][k] = wrap_2pi(X_a_pred[3][k]);
      
     // 6. Copiar ruidos tal cual al estado ampliado predicho
    for (i = 0; i < NW; i++)
        X_a_pred[NX + i][k] = X_a[NX + i][k];

    for (i = 0; i < NV; i++)
        X_a_pred[NX + NW + i][k] = X_a[NX + NW + i][k];
    
    }


   // 7. Media del estado predicho
for (i = 0; i < NX; i++) {
    x_hat_pred[i] = 0.0;
    for (k = 0; k < N; k++)
        x_hat_pred[i] += Wmedia[k] * X_a_pred[i][k];
}
      x_aux[0] = X_a_pred[0][0];  
      x_aux[1] = X_a_pred[1][0];  
      x_aux[2] = X_a_pred[2][0];  
      x_aux[3] = X_a_pred[3][0];  

// 8. Covarianza predicha
for (i = 0; i < NX; i++)
    for (j = 0; j < NX; j++)
        P_pred[i][j] = 0.0;

for (k = 0; k < N; k++) {
    for (i = 0; i < NX; i++)
        dx[i] = X_a_pred[i][k] - x_hat_pred[i];

    for (i = 0; i < NX; i++)
        for (j = 0; j < NX; j++)
            P_pred[i][j] += Wcov[k] * dx[i] * dx[j];
}


for (i = 0; i < NX; i++)
    for (j = 0; j < NX; j++)
        P[i][j] = P_pred[i][j];




//10)Implementamos las mediciones de corrientes más el ruido v
    for (k = 0; k < N; k++)
{
    // Medición de Ia
    Y_a_pred[0][k] =
        X_a_pred[0][k] + X_a_pred[NX + NW + 0][k];

    // Medición de Ib
    Y_a_pred[1][k] =
        X_a_pred[1][k] + X_a_pred[NX + NW + 1][k];
}

    for (i = 0; i < NV; i++) {
    y_hat_pred[i] = 0.0;
    for (k = 0; k < N; k++)
        y_hat_pred[i] += Wmedia[k] * Y_a_pred[i][k];
    }

 //11) covarianza de las mediciones

    for (i = 0; i < NV; i++)
    for (j = 0; j < NV; j++)
        S_cov[i][j] = 0.0;

for (k = 0; k < N; k++)
{
    for (i = 0; i < NV; i++)
        dy[i] = Y_a_pred[i][k] - y_hat_pred[i];

    for (i = 0; i < NV; i++)
        for (j = 0; j < NV; j++)
            S_cov[i][j] += Wcov[k] * dy[i] * dy[j];
}


    //12) covarianzas cruzadas
    for (i = 0; i < NX; i++)
    for (j = 0; j < NV; j++)
        P_xy[i][j] = 0.0;

for (k = 0; k < N; k++)
{
    for (i = 0; i < NX; i++)
        dx[i] = X_a_pred[i][k] - x_hat_pred[i];

    for (i = 0; i < NV; i++)
        dy[i] = Y_a_pred[i][k] - y_hat_pred[i];

    for (i = 0; i < NX; i++)
        for (j = 0; j < NV; j++)
            P_xy[i][j] += Wcov[k] * dx[i] * dy[j];
}


    //13) Ganancias de kalman

    double a = S_cov[0][0];
double b = S_cov[0][1];
double c = S_cov[1][0];
double d = S_cov[1][1];

double det = a*d - b*c;
    if (fabs(det) < 1e-12) det = 1e-12;

double inv_det = 1.0 / det;


S_inv[0][0] =  d * inv_det;
S_inv[0][1] = -b * inv_det;
S_inv[1][0] = -c * inv_det;
S_inv[1][1] =  a * inv_det;

// K = P_xy * S_inv

for (i = 0; i < NX; i++)
    for (j = 0; j < NV; j++)
        K_gain[i][j] = 
            P_xy[i][0] * S_inv[0][j] +
            P_xy[i][1] * S_inv[1][j];

    //14) Actualización del estado

innovation[0] = (*Ia_medido) - y_hat_pred[0];
innovation[1] = (*Ib_medido) - y_hat_pred[1];

for (i = 0; i < NX; i++){
    x_hat[i] =  x_hat_pred[i] + K_gain[i][0] * innovation[0] + K_gain[i][1] * innovation[1];
    x[i] = x_hat[i];
}
   x_hat[3] = wrap_2pi(x_hat[3]);
   x[3] = wrap_2pi(x[3]);


 //15) Actualización de la covarianza
// KS = K * S_cov
for (i = 0; i < NX; i++)
    for (j = 0; j < NV; j++)
        KS[i][j] = K_gain[i][0] * S_cov[0][j] + K_gain[i][1] * S_cov[1][j];

for (i = 0; i < NX; i++)
    for (j = 0; j < NX; j++)
        P[i][j] = P_pred[i][j] -  (KS[i][0] * K_gain[j][0] + KS[i][1] * K_gain[j][1]);
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


