
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
/* ======================= Estado persistente ======================= */
/* ==== PARÁMETROS (ajusta aquí) ==== */


/* ==== ESTADOS ==== */
static double ek2_d=0, ek1_d=0, ek_d=0; /*errores corriente d */
static double ek2_q=0, ek1_q=0, ek_q=0; /*errores corriente q */
static double ek2_w=0, ek1_w=0, ek_w=0; /*errores velocidad */
static double ud=0, uq=0;               /*Salidas voltaje */

/* ==== UTILS ==== */
static double clamp(double v, double vmin, double vmax){
    if(v<vmin) return vmin;
    if(v>vmax) return vmax;
    return v;
}
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
#define u_10_width 1
#define u_11_width 1
#define u_12_width 1
#define u_13_width 1
#define u_14_width 1
#define u_15_width 1
#define y_width 1
#define y_1_width 1
#define y_2_width 1
#define y_3_width 1

/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
 
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Start function
 *
 */
void PI_dq_Start_wrapper(void)
{
/* %%%-SFUNWIZ_wrapper_Start_Changes_BEGIN --- EDIT HERE TO _END */
    ek2_d = ek1_d = ek_d = 0.0;
    ek2_q = ek1_q = ek_q = 0.0;
    ek2_w = ek1_w = ek_w = 0.0;
    ud = 0.0;
    uq = 0.0;
/* %%%-SFUNWIZ_wrapper_Start_Changes_END --- EDIT HERE TO _BEGIN */
}
/*
 * Output function
 *
 */
void PI_dq_Outputs_wrapper(const real_T *Ia,
			const real_T *Ib,
			const real_T *Wm_ref,
			const real_T *Wm_ext,
			const real_T *Kp_current,
			const real_T *Ki_current,
			const real_T *Kd_current,
			const real_T *sample_time_ext,
			const real_T *Vdc_ext,
			const real_T *R_ext,
			const real_T *Ld_ext,
			const real_T *Lq_ext,
			const real_T *Ke_ext,
			const real_T *Theta_ext,
			const real_T *Kp_wm,
			const real_T *Ki_wm,
			real_T *Ud,
			real_T *Uq,
			real_T *Id,
			real_T *Iq)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
/* Parámetros on-line */
    
    
    
    
    
    const double T   = sample_time_ext[0];
    const double Vdc = Vdc_ext[0];
    const double Imax = 1.0; /* 1 A por fase */
    const double Wm_max = 80;

    /* parámetros Kp,Kd y Ki ára corriente y velocidad */
     const double Kp_c  = Kp_current[0];
     const double Ki_c = Ki_current[0];
     const double Kp_w  = Kp_wm[0];
     const double Ki_w = Ki_wm[0];
    /* Parámetros del motor */

    const double R   = R_ext[0];
    const double Ld  = Ld_ext[0];
    const double Lq  = Lq_ext[0];
    const double Ke  = Ke_ext[0];

    /* Lecturas de variables*/
     const double Wm  = Wm_ext[0];
     const double s = sin(Theta_ext[0]);
     const double c = cos(Theta_ext[0]);


      Id[0] =  Ia[0]*c + Ib[0]*s;
      Iq[0] = -Ia[0]*s + Ib[0]*c;
      /* ---- Errores Wm ---- */
    ek2_w = ek1_w;
    ek1_w = ek_w;
        {
        double idr = clamp(Wm_ref[0], -Wm_max, Wm_max);
        ek_w = Wm_ref[0] - Wm;
        }
    /*  Aplicamos errores PID discretos*/
        
     /*K1 y K2 de Wm */
    const double k1_w = (Kp_w*T) + (Ki_w*T*T);
    const double k2_w = -(Kp_w*T);
    const double Iq_ref = (1.0/T)*(k1_w*ek_w + k2_w*ek1_w)/Ke;
    const double Id_ref = 0;


    /* ---- Errores d ---- */
    ek2_d = ek1_d;
    ek1_d = ek_d;
    {
        double idr = clamp(Id_ref, -Imax, Imax);
        ek_d = idr - Id[0];
    }

    /* ---- Errores q ---- */
    ek2_q = ek1_q;
    ek1_q = ek_q;
    {
        double iqr = clamp(Iq_ref, -Imax, Imax);
        ek_q = iqr - Iq[0];
    }

    /* ---- PI incremental (Euler backward) ----
       u_ctrl[k] = u_ctrl[k-1] + (1/T)*(k1*e[k] + k2*e[k-1])
       (PI => Kd=0)
    */
    {   
       
        /*K1 y K2 de Idq */
        const double k1_c = (Kp_c*T) + (Ki_c*T*T); 
        const double k2_c = -(Kp_c*T);

        
        const double du_d = (1.0/T)*(k1_c*ek_d + k2_c*ek1_d);
        const double du_q = (1.0/T)*(k1_c*ek_q + k2_c*ek1_q);

        /* Candidatos del PI (solo parte de control) */
        double ud_ctrl_trial = ud + du_d;
        double uq_ctrl_trial = uq + du_q;

        /* --------- Feed-forward / desacople --------- */
        /* ud_ff = R*id - we*Lq*iq */
        /* uq_ff = R*iq + we*Ld*id + we*Ke */
        {

            const double p = 50;
            const double ud_ff = + p*Wm*Lq*Iq[0];
            const double uq_ff = - p*Wm*Ld*Id[0] - p*Wm*Ke;

            /* Salida total a evaluar/saturar */
            double ud_tot_trial = ud_ctrl_trial + ud_ff;
            double uq_tot_trial = uq_ctrl_trial + uq_ff;

            /* Anti-windup incremental (clamping) decide con salida TOTAL */
            /* Si saturaría y el incremento del PI empuja hacia fuera, no integres */
            if (!((ud_tot_trial >  Vdc && (ud_ctrl_trial - ud) > 0.0) ||
                  (ud_tot_trial < -Vdc && (ud_ctrl_trial - ud) < 0.0)))
            {
                ud = ud_ctrl_trial;
            }

            if (!((uq_tot_trial >  Vdc && (uq_ctrl_trial - uq) > 0.0) ||
                  (uq_tot_trial < -Vdc && (uq_ctrl_trial - uq) < 0.0)))
            {
                uq = uq_ctrl_trial;
            }

            /* Reconstruye salida TOTAL ya con estado actualizado */
            ud_tot_trial = ud + ud_ff;
            uq_tot_trial = uq + uq_ff;

            /* Saturación por eje */
            ud_tot_trial = clamp(ud_tot_trial, -Vdc, Vdc);
            uq_tot_trial = clamp(uq_tot_trial, -Vdc, Vdc);

            /* Saturación vectorial: ||U|| <= Vdc */
            {
                double Umag = sqrt(ud_tot_trial*ud_tot_trial + uq_tot_trial*uq_tot_trial);
                if (Umag > Vdc && Umag > 0.0){
                    double s = Vdc/Umag;
                    ud_tot_trial *= s;
                    uq_tot_trial *= s;
                }
            }

            /* Salidas */
            Ud[0] = ud_tot_trial;
            Uq[0] = uq_tot_trial;
        }
    }
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}


