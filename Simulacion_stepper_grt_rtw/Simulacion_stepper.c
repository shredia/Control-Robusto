/*
 * Simulacion_stepper.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Simulacion_stepper".
 *
 * Model version              : 1.44
 * Simulink Coder version : 24.1 (R2024a) 19-Nov-2023
 * C source code generated on : Mon Nov 24 00:43:23 2025
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Debugging
 * Validation result: Not run
 */

#include "Simulacion_stepper.h"
#include "rtwtypes.h"
#include "Simulacion_stepper_private.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>
#include <emmintrin.h>
#include <float.h>

/* Block signals (default storage) */
B_Simulacion_stepper_T Simulacion_stepper_B;

/* Block states (default storage) */
DW_Simulacion_stepper_T Simulacion_stepper_DW;

/* Real-time model */
static RT_MODEL_Simulacion_stepper_T Simulacion_stepper_M_;
RT_MODEL_Simulacion_stepper_T *const Simulacion_stepper_M =
  &Simulacion_stepper_M_;

/* Forward declaration for local functions */
static real_T Simulacion_stepper_mod(real_T x);
static void rate_scheduler(void);
real_T look1_pbinlxmpw(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T prevIndex[], uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T bpIdx;

  /* Row-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'on'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'on'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    bpIdx = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    uint32_T found;
    uint32_T iLeft;
    uint32_T iRght;

    /* Binary Search using Previous Index */
    bpIdx = prevIndex[0U];
    iLeft = 0U;
    iRght = maxIndex;
    found = 0U;
    while (found == 0U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx - 1U;
        bpIdx = ((bpIdx + iLeft) - 1U) >> 1U;
      } else if (u0 < bp0[bpIdx + 1U]) {
        found = 1U;
      } else {
        iLeft = bpIdx + 1U;
        bpIdx = ((bpIdx + iRght) + 1U) >> 1U;
      }
    }

    frac = (u0 - bp0[bpIdx]) / (bp0[bpIdx + 1U] - bp0[bpIdx]);
  } else {
    bpIdx = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  prevIndex[0U] = bpIdx;

  /* Row-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[bpIdx];
  return (table[bpIdx + 1U] - yL_0d0) * frac + yL_0d0;
}

int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator)
{
  return (((numerator < 0) != (denominator < 0)) && (numerator % denominator !=
           0) ? -1 : 0) + numerator / denominator;
}

/*
 *         This function updates active task flag for each subrate.
 *         The function is called at model base rate, hence the
 *         generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (Simulacion_stepper_M->Timing.TaskCounters.TID[1])++;
  if ((Simulacion_stepper_M->Timing.TaskCounters.TID[1]) > 4) {/* Sample time: [5.0E-6s, 0.0s] */
    Simulacion_stepper_M->Timing.TaskCounters.TID[1] = 0;
  }
}

/*
 * Output and update for atomic system:
 *    '<Root>/MATLAB Function'
 *    '<Root>/MATLAB Function3'
 *    '<Root>/MATLAB Function4'
 */
void Simulacion_stepp_MATLABFunction(real_T rtu_u,
  B_MATLABFunction_Simulacion_s_T *localB)
{
  real_T angle;

  /* :  y = wrap_2pi(u); */
  /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
  if (rtIsNaN(rtu_u + 3.1415926535897931) || rtIsInf(rtu_u + 3.1415926535897931))
  {
    angle = (rtNaN);
  } else if (rtu_u + 3.1415926535897931 == 0.0) {
    angle = 0.0;
  } else {
    boolean_T rEQ0;
    angle = fmod(rtu_u + 3.1415926535897931, 6.2831853071795862);
    rEQ0 = (angle == 0.0);
    if (!rEQ0) {
      real_T q;
      q = fabs((rtu_u + 3.1415926535897931) / 6.2831853071795862);
      rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
    }

    if (rEQ0) {
      angle = 0.0;
    } else if (rtu_u + 3.1415926535897931 < 0.0) {
      angle += 6.2831853071795862;
    }
  }

  /* 'wrap_2pi:3' if angle < 0 */
  if (angle < 0.0) {
    /* 'wrap_2pi:4' angle = angle + 2*pi; */
    angle += 6.2831853071795862;
  }

  /* 'wrap_2pi:7' angle_output = angle - pi; */
  localB->y = angle - 3.1415926535897931;
}

/*
 * Output and update for atomic system:
 *    '<Root>/MATLAB Function1'
 *    '<Root>/MATLAB Function2'
 */
void Simulacion_step_MATLABFunction1(real_T rtu_u,
  B_MATLABFunction1_Simulacion__T *localB)
{
  real_T angle;

  /* :  y = wrap_2pi(u); */
  /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
  if (rtIsNaN(rtu_u + 3.1415926535897931) || rtIsInf(rtu_u + 3.1415926535897931))
  {
    angle = (rtNaN);
  } else if (rtu_u + 3.1415926535897931 == 0.0) {
    angle = 0.0;
  } else {
    boolean_T rEQ0;
    angle = fmod(rtu_u + 3.1415926535897931, 6.2831853071795862);
    rEQ0 = (angle == 0.0);
    if (!rEQ0) {
      real_T q;
      q = fabs((rtu_u + 3.1415926535897931) / 6.2831853071795862);
      rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
    }

    if (rEQ0) {
      angle = 0.0;
    } else if (rtu_u + 3.1415926535897931 < 0.0) {
      angle += 6.2831853071795862;
    }
  }

  /* 'wrap_2pi:3' if angle < 0 */
  if (angle < 0.0) {
    /* 'wrap_2pi:4' angle = angle + 2*pi; */
    angle += 6.2831853071795862;
  }

  /* 'wrap_2pi:7' angle_output = angle - pi; */
  localB->y = angle - 3.1415926535897931;
}

/* Function for MATLAB Function: '<S23>/MATLAB Function' */
static real_T Simulacion_stepper_mod(real_T x)
{
  real_T r;
  if (rtIsNaN(x) || rtIsInf(x)) {
    r = (rtNaN);
  } else if (x == 0.0) {
    r = 0.0;
  } else {
    boolean_T rEQ0;
    r = fmod(x, 6.2831853071795862);
    rEQ0 = (r == 0.0);
    if (!rEQ0) {
      real_T q;
      q = fabs(x / 6.2831853071795862);
      rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
    }

    if (rEQ0) {
      r = 0.0;
    } else if (x < 0.0) {
      r += 6.2831853071795862;
    }
  }

  return r;
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (rtIsInf(u1)) {
    y = u0;
  } else if ((u1 != 0.0) && (u1 != trunc(u1))) {
    real_T q;
    q = fabs(u0 / u1);
    if (!(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q)) {
      y = 0.0 * u0;
    } else {
      y = fmod(u0, u1);
    }
  } else {
    y = fmod(u0, u1);
  }

  return y;
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T y;
  y = u0;
  if (u1 == 0.0) {
    if (u0 == 0.0) {
      y = u1;
    }
  } else if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (u0 == 0.0) {
    y = 0.0 / u1;
  } else if (rtIsInf(u1)) {
    if ((u1 < 0.0) != (u0 < 0.0)) {
      y = u1;
    }
  } else {
    boolean_T yEq;
    y = fmod(u0, u1);
    yEq = (y == 0.0);
    if ((!yEq) && (u1 > floor(u1))) {
      real_T q;
      q = fabs(u0 / u1);
      yEq = !(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q);
    }

    if (yEq) {
      y = u1 * 0.0;
    } else if ((u0 < 0.0) != (u1 < 0.0)) {
      y += u1;
    }
  }

  return y;
}

/* Model step function */
void Simulacion_stepper_step(void)
{
  __m128d tmp_b;
  __m128d tmp_d;
  real_T b[441];
  real_T b_b[441];
  real_T X[210];
  real_T S_0[100];
  real_T X_minus[84];
  real_T Xs[84];
  real_T Y_minus[42];
  real_T Ys[42];
  real_T Ys_0[42];
  real_T A_hat_0[25];
  real_T A_hat_2[25];
  real_T P_pred_0[25];
  real_T Qk_0[25];
  real_T A_hat[16];
  real_T A_hat_1[16];
  real_T P_pred[16];
  real_T Qk[16];
  real_T Lk_0[10];
  real_T a[10];
  real_T B[8];
  real_T Lk[8];
  real_T tmp_0[5];
  real_T tmp_1[5];
  real_T tmp_2[5];
  real_T x_pred_0[5];
  real_T Rk[4];
  real_T Sk[4];
  real_T x_pred[4];
  real_T tmp[2];
  real_T tmp_c[2];
  real_T y_pred[2];
  real_T RL;
  real_T S;
  real_T Th_e;
  real_T dos;
  real_T tmp_3;
  real_T tmp_4;
  real_T tmp_5;
  real_T tmp_6;
  real_T tmp_7;
  real_T tmp_8;
  real_T tmp_9;
  real_T tmp_a;
  real_T uno;
  int32_T Lk_tmp;
  int32_T b_ix;
  int32_T b_iy;
  int32_T e;
  int32_T ia;
  int32_T idxAjj;
  int32_T iy;
  int32_T jm1;
  int32_T r1;
  int32_T r2;
  int8_T a_0;
  int8_T a_1;
  int8_T b_b_0;
  int8_T b_b_1;
  static const int8_T a_2[8] = { 1, 0, 0, 1, 0, 0, 0, 0 };

  static const int8_T b_b_2[8] = { 1, 0, 0, 0, 0, 1, 0, 0 };

  static const real_T b_0[5] = { 0.01, 0.01, 1.0E-6, 1.0E-10, 1.0E-10 };

  static const int8_T c[5] = { 0, 0, 1, 0, 0 };

  static const int8_T a_3[10] = { 1, 0, 0, 1, 0, 0, 0, 0, 0, 0 };

  static const int8_T b_b_3[10] = { 1, 0, 0, 0, 0, 0, 1, 0, 0, 0 };

  boolean_T exitg1;
  if (Simulacion_stepper_M->Timing.TaskCounters.TID[1] == 0) {
    /* DiscreteIntegrator: '<S90>/Discrete-Time Integrator' */
    Simulacion_stepper_B.IphA[0] =
      Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0];
    Simulacion_stepper_B.IphA[1] =
      Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[1];

    /* S-Function (sfun_spssw_discc): '<S97>/State-Space' incorporates:
     *  Constant: '<S29>/DC'
     *  Constant: '<S30>/DC'
     *  Constant: '<S99>/SwitchCurrents'
     */

    /* S-Function block: <S97>/State-Space */
    {
      real_T accum;

      /* Circuit has switches */
      int_T *switch_status = (int_T*)
        Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS;
      int_T *switch_status_init = (int_T*)
        Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;
      int_T *SwitchChange = (int_T*)
        Simulacion_stepper_DW.StateSpace_PWORK.SW_CHG;
      int_T *gState = (int_T*) Simulacion_stepper_DW.StateSpace_PWORK.G_STATE;
      real_T *yswitch = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.Y_SWITCH;
      int_T *switchTypes = (int_T*)
        Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_TYPES;
      int_T *idxOutSw = (int_T*)
        Simulacion_stepper_DW.StateSpace_PWORK.IDX_OUT_SW;
      real_T *DxCol = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.DX_COL;
      real_T *tmp2 = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.TMP2;
      real_T *uswlast = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.USWLAST;
      int_T newState;
      int_T swChanged = 0;
      int loopsToDo = 20;
      real_T temp;

      /* keep an initial copy of switch_status*/
      memcpy(switch_status_init, switch_status, 16 * sizeof(int_T));
      memcpy(uswlast, &Simulacion_stepper_B.StateSpace_o1[0], 16*sizeof(real_T));
      do {
        if (loopsToDo == 1) {          /* Need to reset some variables: */
          swChanged = 0;

          /* return to the original switch status*/
          {
            int_T i1;
            for (i1=0; i1 < 16; i1++) {
              swChanged = ((SwitchChange[i1] = switch_status_init[i1] -
                            switch_status[i1]) != 0) ? 1 : swChanged;
              switch_status[i1] = switch_status_init[i1];
            }
          }
        } else {
          /*
           * Compute outputs:
           * ---------------
           */
          real_T *Ds = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.DS;

          {
            int_T i1;
            real_T *y0 = &Simulacion_stepper_B.StateSpace_o1[0];
            for (i1=0; i1 < 22; i1++) {
              accum = 0.0;

              {
                int_T i2;
                const real_T *u0 = &Simulacion_stepper_P.SwitchCurrents_Value[0];
                for (i2=0; i2 < 16; i2++) {
                  accum += *(Ds++) * u0[i2];
                }

                accum += *(Ds++) * Simulacion_stepper_B.IphA[0];
                accum += *(Ds++) * Simulacion_stepper_B.IphA[1];
                accum += *(Ds++) * Simulacion_stepper_P.Vdc;
                accum += *(Ds++) * Simulacion_stepper_P.Vdc;
              }

              y0[i1] = accum;
            }
          }

          swChanged = 0;

          {
            int_T i1;
            real_T *y0 = &Simulacion_stepper_B.StateSpace_o1[0];
            for (i1=0; i1 < 16; i1++) {
              switch (switchTypes[i1]) {
               case 1 :                /* Ideal switch */
                newState = gState[i1] > 0 ? 1 : 0;
                break;

               case 3 :                /* Diodes */
                newState = y0[i1] > 0.0 ? 1 : ((y0[i1] < 0.0) ? 0 :
                  switch_status[i1]);
                break;
              }

              swChanged = ((SwitchChange[i1] = newState - switch_status[i1]) !=
                           0) ? 1 : swChanged;
              switch_status[i1] = newState;/* Keep new state */
            }
          }
        }

        /*
         * Compute new As, Bs, Cs and Ds matrixes:
         * --------------------------------------
         */
        if (swChanged) {
          real_T *Ds = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.DS;
          real_T a1;

          {
            int_T i1;
            for (i1=0; i1 < 16; i1++) {
              if (SwitchChange[i1] != 0) {
                a1 = yswitch[i1]*SwitchChange[i1];
                temp = 1/(1-Ds[i1*21]*a1);

                {
                  int_T i2;
                  for (i2=0; i2 < 22; i2++) {
                    DxCol[i2]= Ds[i2 * 20 + i1]*temp*a1;
                  }
                }

                DxCol[i1] = temp;

                /* Copy row nSw of Ds into tmp2 and zero it out in Ds */
                memcpy(tmp2, &Ds[i1 * 20], 20 * sizeof(real_T));
                memset(&Ds[i1 * 20], '\0', 20 * sizeof(real_T));

                /* Cs = Cs + DxCol * tmp1, Ds = Ds + DxCol * tmp2 *******************/
                {
                  int_T i2;
                  for (i2=0; i2 < 22; i2++) {
                    a1 = DxCol[i2];

                    {
                      int_T i3;
                      for (i3=0; i3 < 20; i3++) {
                        Ds[i2 * 20 + i3] += a1 * tmp2[i3];
                      }
                    }
                  }
                }
              }
            }
          }
        }                              /* if (swChanged) */
      } while (swChanged > 0 && --loopsToDo > 0);

      if (loopsToDo == 0) {
        real_T *Ds = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.DS;

        {
          int_T i1;
          real_T *y0 = &Simulacion_stepper_B.StateSpace_o1[0];
          for (i1=0; i1 < 22; i1++) {
            accum = 0.0;

            {
              int_T i2;
              const real_T *u0 = &Simulacion_stepper_P.SwitchCurrents_Value[0];
              for (i2=0; i2 < 16; i2++) {
                accum += *(Ds++) * u0[i2];
              }

              accum += *(Ds++) * Simulacion_stepper_B.IphA[0];
              accum += *(Ds++) * Simulacion_stepper_B.IphA[1];
              accum += *(Ds++) * Simulacion_stepper_P.Vdc;
              accum += *(Ds++) * Simulacion_stepper_P.Vdc;
            }

            y0[i1] = accum;
          }
        }
      }

      /* Output new switches states */
      {
        int_T i1;
        real_T *y1 = &Simulacion_stepper_B.StateSpace_o2[0];
        for (i1=0; i1 < 16; i1++) {
          y1[i1] = (real_T)switch_status[i1];
        }
      }
    }

    /* Gain: '<S87>/eliminate warning with bus selector' */
    tmp_d = _mm_mul_pd(_mm_set1_pd
                       (Simulacion_stepper_P.eliminatewarningwithbusselector),
                       _mm_loadu_pd(&Simulacion_stepper_B.StateSpace_o1[16]));

    /* Gain: '<S87>/eliminate warning with bus selector' */
    _mm_storeu_pd(&Simulacion_stepper_B.VphV[0], tmp_d);

    /* Gain: '<S2>/do not delete this gain' incorporates:
     *  Gain: '<S1>/do not delete this gain'
     */
    tmp_d = _mm_mul_pd(_mm_set_pd
                       (Simulacion_stepper_P.donotdeletethisgain_Gain_l,
                        Simulacion_stepper_P.donotdeletethisgain_Gain),
                       _mm_loadu_pd(&Simulacion_stepper_B.StateSpace_o1[20]));
    _mm_storeu_pd(&tmp_c[0], tmp_d);

    /* Gain: '<S1>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain = tmp_c[0];

    /* Gain: '<S2>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain_d = tmp_c[1];

    /* DiscreteIntegrator: '<S89>/Discrete-Time Integrator1' */
    Simulacion_stepper_B.DiscreteTimeIntegrator1 =
      Simulacion_stepper_DW.DiscreteTimeIntegrator1_DSTATE;

    /* DiscreteIntegrator: '<S89>/Discrete-Time Integrator' */
    Simulacion_stepper_B.DiscreteTimeIntegrator =
      Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE_o;

    /* MATLAB Function: '<Root>/MATLAB Function4' */
    Simulacion_stepp_MATLABFunction(Simulacion_stepper_B.DiscreteTimeIntegrator,
      &Simulacion_stepper_B.sf_MATLABFunction4);

    /* Gain: '<Root>/Gain2' */
    Simulacion_stepper_B.Gain2 = Simulacion_stepper_P.P *
      Simulacion_stepper_B.sf_MATLABFunction4.y;

    /* Gain: '<Root>/Gain5' */
    Simulacion_stepper_B.Gain5 = Simulacion_stepper_P.Gain5_Gain *
      Simulacion_stepper_B.Gain2;

    /* Trigonometry: '<Root>/Sin' */
    Simulacion_stepper_B.Sin = sin(Simulacion_stepper_B.Gain5);

    /* Gain: '<Root>/Gain4' */
    Simulacion_stepper_B.Gain4 = Simulacion_stepper_P.Tdm *
      Simulacion_stepper_B.Sin;

    /* Sum: '<Root>/Sum3' incorporates:
     *  Constant: '<Root>/Constant'
     */
    Simulacion_stepper_B.Sum3 = Simulacion_stepper_P.Constant_Value +
      Simulacion_stepper_B.Gain4;

    /* S-Function (fcgen): '<Root>/Function-Call Generator1' incorporates:
     *  SubSystem: '<Root>/Triggered Subsystem'
     */
    /* S-Function (PI_dq): '<S22>/C//C++ Code Block' incorporates:
     *  Constant: '<Root>/Constant8'
     *  Constant: '<S22>/Constant'
     *  Constant: '<S22>/Constant1'
     *  Constant: '<S22>/Constant10'
     *  Constant: '<S22>/Constant11'
     *  Constant: '<S22>/Constant2'
     *  Constant: '<S22>/Constant3'
     *  Constant: '<S22>/Constant4'
     *  Constant: '<S22>/Constant5'
     *  Constant: '<S22>/Constant6'
     *  Constant: '<S22>/Constant7'
     *  Constant: '<S22>/Constant8'
     *  Constant: '<S22>/Constant9'
     */
    S = -Simulacion_stepper_P.Kt;
    PI_dq_Outputs_wrapper(&Simulacion_stepper_B.donotdeletethisgain,
                          &Simulacion_stepper_B.donotdeletethisgain_d,
                          &Simulacion_stepper_P.Constant8_Value,
                          &Simulacion_stepper_B.DiscreteTimeIntegrator1,
                          &Simulacion_stepper_P.Kp_corriente,
                          &Simulacion_stepper_P.Ki_corriente,
                          &Simulacion_stepper_P.Constant4_Value,
                          &Simulacion_stepper_P.sample_time,
                          &Simulacion_stepper_P.Vdc, &Simulacion_stepper_P.R,
                          &Simulacion_stepper_P.L, &Simulacion_stepper_P.L, &S,
                          &Simulacion_stepper_B.Gain2,
                          &Simulacion_stepper_P.Kp_w, &Simulacion_stepper_P.Ki_w,
                          &Simulacion_stepper_B.Sum3, &Simulacion_stepper_P.B,
                          &Simulacion_stepper_B.CCCodeBlock_o1,
                          &Simulacion_stepper_B.CCCodeBlock_o2,
                          &Simulacion_stepper_B.CCCodeBlock_o3,
                          &Simulacion_stepper_B.CCCodeBlock_o4,
                          &Simulacion_stepper_B.CCCodeBlock_o5);

    /* End of Outputs for S-Function (fcgen): '<Root>/Function-Call Generator1' */
    /* Fcn: '<S18>/alpha' */
    Simulacion_stepper_B.alpha = cos(0.0 - Simulacion_stepper_B.Gain2) *
      Simulacion_stepper_B.CCCodeBlock_o1 + sin(0.0 - Simulacion_stepper_B.Gain2)
      * Simulacion_stepper_B.CCCodeBlock_o2;

    /* Fcn: '<S18>/beta' */
    Simulacion_stepper_B.beta = -sin(0.0 - Simulacion_stepper_B.Gain2) *
      Simulacion_stepper_B.CCCodeBlock_o1 + cos(0.0 - Simulacion_stepper_B.Gain2)
      * Simulacion_stepper_B.CCCodeBlock_o2;

    /* S-Function (fcgen): '<Root>/Function-Call Generator' incorporates:
     *  SubSystem: '<Root>/Triggered Subsystem1'
     */
    /* SignalConversion generated from: '<S92>/ SFunction ' incorporates:
     *  MATLAB Function: '<S23>/MATLAB Function'
     */
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_f[0] =
      Simulacion_stepper_B.alpha;
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_f[1] =
      Simulacion_stepper_B.beta;

    /* SignalConversion generated from: '<S92>/ SFunction ' incorporates:
     *  MATLAB Function: '<S23>/MATLAB Function'
     */
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_o[0] =
      Simulacion_stepper_B.donotdeletethisgain;
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_o[1] =
      Simulacion_stepper_B.donotdeletethisgain_d;

    /* MATLAB Function: '<S23>/MATLAB Function' incorporates:
     *  Constant: '<S23>/Constant'
     *  Constant: '<S23>/Constant19'
     *  Constant: '<S23>/Constant20'
     *  Constant: '<S23>/Constant21'
     *  Constant: '<S23>/Constant22'
     *  Constant: '<S23>/Constant23'
     *  Constant: '<S23>/Constant24'
     */
    /* :  if isempty(x_hat) */
    /* :  Qk = diag([1e-2 1e-2 1e-6 1e-8]); */
    memset(&Qk[0], 0, sizeof(real_T) << 4U);

    /* :  Rk = diag([5e-2 5e-2]); */
    Qk[0] = 0.01;
    Qk[5] = 0.01;
    Rk[1] = 0.0;
    Qk[10] = 1.0E-6;
    Rk[2] = 0.0;
    Qk[15] = 1.0E-8;
    Rk[0] = 0.05;
    Rk[3] = 0.05;

    /* :  Ia = x_hat(1); */
    /* :  Ib = x_hat(2); */
    /* :  Wm = x_hat(3); */
    /* :  Th_m = x_hat(4); */
    /* :  Th_e = Nr*Th_m; */
    Th_e = Simulacion_stepper_P.P * Simulacion_stepper_DW.x_hat_l[3];

    /* :  S = sin(Th_e); */
    S = sin(Th_e);

    /* :  C = cos(Th_e); */
    Th_e = cos(Th_e);

    /* :  Va = U(1); */
    /* :  Vb = U(2); */
    /* :  dx = zeros(4,1); */
    /* :  dx(1) = (Va - R*Ia + Kt*Wm*S)/L; */
    x_pred[0] = ((Simulacion_stepper_B.TmpSignalConversionAtSFunctio_f[0] -
                  Simulacion_stepper_P.R * Simulacion_stepper_DW.x_hat_l[0]) +
                 Simulacion_stepper_P.Kt * Simulacion_stepper_DW.x_hat_l[2] * S)
      / Simulacion_stepper_P.L;

    /* :  dx(2) = (Vb - R*Ib - Kt*Wm*C)/L; */
    x_pred[1] = ((Simulacion_stepper_B.TmpSignalConversionAtSFunctio_f[1] -
                  Simulacion_stepper_P.R * Simulacion_stepper_DW.x_hat_l[1]) -
                 Simulacion_stepper_P.Kt * Simulacion_stepper_DW.x_hat_l[2] *
                 Th_e) / Simulacion_stepper_P.L;

    /* :  dx(3) = (Kt*(Ib*C-Ia*S)-B*Wm)/J; */
    x_pred[2] = ((Simulacion_stepper_DW.x_hat_l[1] * Th_e -
                  Simulacion_stepper_DW.x_hat_l[0] * S) *
                 Simulacion_stepper_P.Kt - Simulacion_stepper_P.B *
                 Simulacion_stepper_DW.x_hat_l[2]) / Simulacion_stepper_P.J;

    /* :  dx(4) = Wm; */
    x_pred[3] = Simulacion_stepper_DW.x_hat_l[2];

    /* :  x_pred = x_hat + Ts*dx; */
    tmp_d = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(Simulacion_stepper_P.sample_time),
      _mm_loadu_pd(&x_pred[0])), _mm_loadu_pd(&Simulacion_stepper_DW.x_hat_l[0]));
    _mm_storeu_pd(&x_pred[0], tmp_d);
    tmp_d = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(Simulacion_stepper_P.sample_time),
      _mm_loadu_pd(&x_pred[2])), _mm_loadu_pd(&Simulacion_stepper_DW.x_hat_l[2]));
    _mm_storeu_pd(&x_pred[2], tmp_d);

    /* :  x_pred(4) = wrap_2pi(x_pred(4)); */
    /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
    RL = Simulacion_stepper_mod(x_pred[3] + 3.1415926535897931);

    /* 'wrap_2pi:3' if angle < 0 */
    if (RL < 0.0) {
      /* 'wrap_2pi:4' angle = angle + 2*pi; */
      RL += 6.2831853071795862;
    }

    /* 'wrap_2pi:7' angle_output = angle - pi; */
    x_pred[3] = RL - 3.1415926535897931;

    /* :  RL = R/L; */
    RL = Simulacion_stepper_P.R / Simulacion_stepper_P.L;

    /* :  uno = Kt*S; */
    uno = Simulacion_stepper_P.Kt * S;

    /* :  dos = Kt*C; */
    dos = Simulacion_stepper_P.Kt * Th_e;

    /* :  A = [-RL 0   (uno/L) (dos*Nr*Wm/L); */
    /* :   0  -RL   -dos/L    uno*Nr*Wm/L; */
    /* :   -uno/J dos/J -B/J -Kt*Nr*(Ia*C+Ib*S)/J; */
    /* :   0 0 1 0]; */
    /* :  A_hat = eye(4) + Ts*A; */
    memset(&P_pred[0], 0, sizeof(real_T) << 4U);
    P_pred[0] = 1.0;
    P_pred[5] = 1.0;
    P_pred[10] = 1.0;
    P_pred[15] = 1.0;
    tmp_3 = Simulacion_stepper_P.sample_time * -RL;
    tmp_4 = Simulacion_stepper_P.sample_time * 0.0;
    tmp_5 = uno / Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    tmp_6 = dos * Simulacion_stepper_P.P * Simulacion_stepper_DW.x_hat_l[2] /
      Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    tmp_7 = Simulacion_stepper_P.sample_time * 0.0;
    RL = Simulacion_stepper_P.sample_time * -RL;
    tmp_8 = -dos / Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    tmp_9 = uno * Simulacion_stepper_P.P * Simulacion_stepper_DW.x_hat_l[2] /
      Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    uno = -uno / Simulacion_stepper_P.J * Simulacion_stepper_P.sample_time;
    dos = dos / Simulacion_stepper_P.J * Simulacion_stepper_P.sample_time;
    tmp_a = -Simulacion_stepper_P.B / Simulacion_stepper_P.J *
      Simulacion_stepper_P.sample_time;
    S = (Simulacion_stepper_DW.x_hat_l[0] * Th_e +
         Simulacion_stepper_DW.x_hat_l[1] * S) * (-Simulacion_stepper_P.Kt *
      Simulacion_stepper_P.P) / Simulacion_stepper_P.J *
      Simulacion_stepper_P.sample_time;
    tmp_d = _mm_add_pd(_mm_set_pd(tmp_7, P_pred[0]), _mm_set_pd(P_pred[1], tmp_3));
    _mm_storeu_pd(&A_hat[0], tmp_d);
    A_hat[2] = uno + P_pred[2];
    A_hat[3] = Simulacion_stepper_P.sample_time * 0.0 + P_pred[3];
    tmp_d = _mm_add_pd(_mm_set_pd(RL, P_pred[4]), _mm_set_pd(P_pred[5], tmp_4));
    _mm_storeu_pd(&A_hat[4], tmp_d);
    A_hat[6] = dos + P_pred[6];
    A_hat[7] = Simulacion_stepper_P.sample_time * 0.0 + P_pred[7];
    tmp_d = _mm_add_pd(_mm_set_pd(tmp_8, P_pred[8]), _mm_set_pd(P_pred[9], tmp_5));
    _mm_storeu_pd(&A_hat[8], tmp_d);
    tmp_d = _mm_add_pd(_mm_set_pd(P_pred[11], tmp_a), _mm_set_pd
                       (Simulacion_stepper_P.sample_time, P_pred[10]));
    _mm_storeu_pd(&A_hat[10], tmp_d);
    tmp_d = _mm_add_pd(_mm_set_pd(tmp_9, P_pred[12]), _mm_set_pd(P_pred[13],
      tmp_6));
    _mm_storeu_pd(&A_hat[12], tmp_d);
    A_hat[14] = S + P_pred[14];
    A_hat[15] = Simulacion_stepper_P.sample_time * 0.0 + P_pred[15];

    /* :  A_hat_T = transpose(A_hat); */
    /* :  P_pred = A_hat * P * A_hat_T + Qk; */
    for (r2 = 0; r2 < 4; r2++) {
      Th_e = A_hat[r2];
      tmp_3 = A_hat[r2 + 4];
      tmp_4 = A_hat[r2 + 8];
      tmp_5 = A_hat[r2 + 12];
      for (r1 = 0; r1 < 4; r1++) {
        jm1 = r1 << 2;
        S = Simulacion_stepper_DW.P_n[jm1] * Th_e;
        S += Simulacion_stepper_DW.P_n[jm1 + 1] * tmp_3;
        S += Simulacion_stepper_DW.P_n[jm1 + 2] * tmp_4;
        S += Simulacion_stepper_DW.P_n[jm1 + 3] * tmp_5;
        A_hat_1[r2 + jm1] = S;
      }

      Th_e = A_hat_1[r2];
      tmp_3 = A_hat_1[r2 + 4];
      tmp_4 = A_hat_1[r2 + 8];
      tmp_5 = A_hat_1[r2 + 12];
      for (r1 = 0; r1 < 4; r1++) {
        S = Th_e * A_hat[r1];
        S += A_hat[r1 + 4] * tmp_3;
        S += A_hat[r1 + 8] * tmp_4;
        S += A_hat[r1 + 12] * tmp_5;
        jm1 = (r1 << 2) + r2;
        P_pred[jm1] = Qk[jm1] + S;
      }
    }

    /* :  Ck = [1 0 0 0; */
    /* :        0 1 0 0]; */
    /* :  Ck_T = transpose(Ck); */
    /* :  y_pred = Ck*x_pred; */
    /* :  Sk = Ck*P_pred*Ck_T + Rk; */
    /* :  Lk = P_pred * Ck_T /Sk; */
    for (r2 = 0; r2 < 2; r2++) {
      a_0 = a_2[r2];
      a_1 = a_2[r2 + 2];
      for (r1 = 0; r1 < 4; r1++) {
        jm1 = r1 << 2;
        Th_e = P_pred[jm1] * (real_T)a_0;
        Th_e += P_pred[jm1 + 1] * (real_T)a_1;
        Th_e += P_pred[jm1 + 2] * 0.0;
        Th_e += P_pred[jm1 + 3] * 0.0;
        Lk[r2 + (r1 << 1)] = Th_e;
      }

      Th_e = Lk[r2];
      tmp_3 = Lk[r2 + 2];
      tmp_4 = Lk[r2 + 4];
      tmp_5 = Lk[r2 + 6];
      for (r1 = 0; r1 < 2; r1++) {
        jm1 = r1 << 2;
        S = (real_T)b_b_2[jm1] * Th_e;
        S += (real_T)b_b_2[jm1 + 1] * tmp_3;
        S += (real_T)b_b_2[jm1 + 2] * tmp_4;
        S += (real_T)b_b_2[jm1 + 3] * tmp_5;
        jm1 = (r1 << 1) + r2;
        Sk[jm1] = Rk[jm1] + S;
      }

      jm1 = r2 << 2;
      a_0 = b_b_2[jm1];
      a_1 = b_b_2[jm1 + 1];
      b_b_0 = b_b_2[jm1 + 2];
      b_b_1 = b_b_2[jm1 + 3];
      for (r1 = 0; r1 <= 2; r1 += 2) {
        tmp_d = _mm_loadu_pd(&P_pred[r1]);
        tmp_d = _mm_mul_pd(_mm_set1_pd(a_0), tmp_d);
        tmp_b = _mm_loadu_pd(&P_pred[r1 + 4]);
        tmp_b = _mm_mul_pd(_mm_set1_pd(a_1), tmp_b);
        tmp_d = _mm_add_pd(tmp_b, tmp_d);
        tmp_b = _mm_loadu_pd(&P_pred[r1 + 8]);
        tmp_b = _mm_mul_pd(_mm_set1_pd(b_b_0), tmp_b);
        tmp_d = _mm_add_pd(tmp_b, tmp_d);
        tmp_b = _mm_loadu_pd(&P_pred[r1 + 12]);
        tmp_b = _mm_mul_pd(_mm_set1_pd(b_b_1), tmp_b);
        tmp_d = _mm_add_pd(tmp_b, tmp_d);
        _mm_storeu_pd(&B[r1 + jm1], tmp_d);
      }
    }

    if (fabs(Sk[1]) > fabs(Sk[0])) {
      r1 = 1;
      r2 = 0;
    } else {
      r1 = 0;
      r2 = 1;
    }

    S = Sk[r2] / Sk[r1];
    tmp_3 = Sk[r1 + 2];
    Th_e = Sk[r2 + 2] - tmp_3 * S;
    Lk_tmp = r1 << 2;
    Lk[Lk_tmp] = B[0] / Sk[r1];
    idxAjj = r2 << 2;
    Lk[idxAjj] = (B[4] - Lk[Lk_tmp] * tmp_3) / Th_e;
    Lk[Lk_tmp] -= Lk[idxAjj] * S;
    r2 = Lk_tmp + 1;
    Lk[r2] = B[1] / Sk[r1];
    jm1 = idxAjj + 1;
    Lk[jm1] = (B[5] - Lk[r2] * tmp_3) / Th_e;
    Lk[r2] -= Lk[jm1] * S;
    r2 = Lk_tmp + 2;
    Lk[r2] = B[2] / Sk[r1];
    jm1 = idxAjj + 2;
    Lk[jm1] = (B[6] - Lk[r2] * tmp_3) / Th_e;
    Lk[r2] -= Lk[jm1] * S;
    Lk_tmp += 3;
    Lk[Lk_tmp] = B[3] / Sk[r1];
    idxAjj += 3;
    Lk[idxAjj] = (B[7] - Lk[Lk_tmp] * tmp_3) / Th_e;
    Lk[Lk_tmp] -= Lk[idxAjj] * S;

    /* :  y_meas = [X(1);X(2)]; */
    /* :  innov = y_meas - y_pred; */
    /* :  x_corr = x_pred + Lk*(innov); */
    y_pred[0] = Simulacion_stepper_B.TmpSignalConversionAtSFunctio_o[0];
    y_pred[1] = Simulacion_stepper_B.TmpSignalConversionAtSFunctio_o[1];
    tmp_3 = x_pred[0];
    S = x_pred[1];
    tmp_4 = x_pred[2];
    tmp_5 = x_pred[3];
    for (r2 = 0; r2 < 2; r2++) {
      Th_e = (real_T)a_2[r2] * tmp_3;
      Th_e += (real_T)a_2[r2 + 2] * S;
      Th_e += 0.0 * tmp_4;
      Th_e += 0.0 * tmp_5;
      tmp[r2] = y_pred[r2] - Th_e;
    }

    S = tmp[0];
    Th_e = tmp[1];
    for (r2 = 0; r2 <= 2; r2 += 2) {
      tmp_d = _mm_loadu_pd(&Lk[r2]);
      tmp_d = _mm_mul_pd(tmp_d, _mm_set1_pd(S));
      tmp_b = _mm_loadu_pd(&Lk[r2 + 4]);
      tmp_b = _mm_mul_pd(tmp_b, _mm_set1_pd(Th_e));
      tmp_d = _mm_add_pd(tmp_b, tmp_d);
      tmp_b = _mm_loadu_pd(&x_pred[r2]);
      tmp_d = _mm_add_pd(tmp_b, tmp_d);
      _mm_storeu_pd(&Simulacion_stepper_B.x_corr_a[r2], tmp_d);
    }

    /* :  x_corr(4) = wrap_2pi(x_corr(4)); */
    /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
    S = Simulacion_stepper_mod(Simulacion_stepper_B.x_corr_a[3] +
      3.1415926535897931);

    /* 'wrap_2pi:3' if angle < 0 */
    if (S < 0.0) {
      /* 'wrap_2pi:4' angle = angle + 2*pi; */
      S += 6.2831853071795862;
    }

    /* 'wrap_2pi:7' angle_output = angle - pi; */
    Simulacion_stepper_B.x_corr_a[3] = S - 3.1415926535897931;

    /* :  x_hat = x_corr; */
    /* :  Lk_T = transpose(Lk); */
    /* :  P = P_pred - Lk*Sk*Lk_T; */
    for (r1 = 0; r1 < 4; r1++) {
      Simulacion_stepper_DW.x_hat_l[r1] = Simulacion_stepper_B.x_corr_a[r1];
      S = Lk[r1];
      Th_e = S * Sk[0];
      tmp_4 = Lk[r1 + 4];
      Th_e += tmp_4 * Sk[1];
      tmp_3 = Th_e;
      Th_e = S * Sk[2];
      Th_e += tmp_4 * Sk[3];
      for (r2 = 0; r2 < 4; r2++) {
        S = tmp_3 * Lk[r2];
        S += Lk[r2 + 4] * Th_e;
        jm1 = (r2 << 2) + r1;
        Simulacion_stepper_DW.P_n[jm1] = P_pred[jm1] - S;
      }
    }

    /* :  P_T = transpose(P); */
    /* :  P = 0.5 * (P + P_T); */
    for (r2 = 0; r2 <= 2; r2 += 2) {
      r1 = r2 << 2;
      jm1 = (r2 + 1) << 2;
      tmp_d = _mm_set1_pd(0.5);
      tmp_b = _mm_mul_pd(_mm_add_pd(_mm_set_pd(Simulacion_stepper_DW.P_n[jm1],
        Simulacion_stepper_DW.P_n[r1]), _mm_loadu_pd
        (&Simulacion_stepper_DW.P_n[r2])), tmp_d);
      _mm_storeu_pd(&tmp_c[0], tmp_b);
      Qk[r1] = tmp_c[0];
      Qk[jm1] = tmp_c[1];
      Lk_tmp = r1 + 1;
      idxAjj = jm1 + 1;
      tmp_b = _mm_mul_pd(_mm_add_pd(_mm_set_pd(Simulacion_stepper_DW.P_n[idxAjj],
        Simulacion_stepper_DW.P_n[Lk_tmp]), _mm_loadu_pd
        (&Simulacion_stepper_DW.P_n[r2 + 4])), tmp_d);
      _mm_storeu_pd(&tmp_c[0], tmp_b);
      Qk[Lk_tmp] = tmp_c[0];
      Qk[idxAjj] = tmp_c[1];
      Lk_tmp = r1 + 2;
      idxAjj = jm1 + 2;
      tmp_b = _mm_mul_pd(_mm_add_pd(_mm_set_pd(Simulacion_stepper_DW.P_n[idxAjj],
        Simulacion_stepper_DW.P_n[Lk_tmp]), _mm_loadu_pd
        (&Simulacion_stepper_DW.P_n[r2 + 8])), tmp_d);
      _mm_storeu_pd(&tmp_c[0], tmp_b);
      Qk[Lk_tmp] = tmp_c[0];
      Qk[idxAjj] = tmp_c[1];
      r1 += 3;
      jm1 += 3;
      tmp_d = _mm_mul_pd(_mm_add_pd(_mm_set_pd(Simulacion_stepper_DW.P_n[jm1],
        Simulacion_stepper_DW.P_n[r1]), _mm_loadu_pd
        (&Simulacion_stepper_DW.P_n[r2 + 12])), tmp_d);
      _mm_storeu_pd(&tmp_c[0], tmp_d);
      Qk[r1] = tmp_c[0];
      Qk[jm1] = tmp_c[1];
    }

    memcpy(&Simulacion_stepper_DW.P_n[0], &Qk[0], sizeof(real_T) << 4U);

    /* SignalConversion generated from: '<S93>/ SFunction ' incorporates:
     *  MATLAB Function: '<S23>/MATLAB Function1'
     */
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_n[0] =
      Simulacion_stepper_B.alpha;
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_n[1] =
      Simulacion_stepper_B.beta;

    /* SignalConversion generated from: '<S93>/ SFunction ' incorporates:
     *  MATLAB Function: '<S23>/MATLAB Function1'
     */
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_h[0] =
      Simulacion_stepper_B.donotdeletethisgain;
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_h[1] =
      Simulacion_stepper_B.donotdeletethisgain_d;

    /* MATLAB Function: '<S23>/MATLAB Function1' incorporates:
     *  Constant: '<S23>/Constant25'
     *  Constant: '<S23>/Constant26'
     *  Constant: '<S23>/Constant27'
     *  Constant: '<S23>/Constant28'
     *  Constant: '<S23>/Constant29'
     *  Constant: '<S23>/Constant30'
     *  Constant: '<S23>/Constant31'
     */
    /* :  if isempty(x_hat) */
    /* :  Qk = diag([1e-2 1e-2 1e-6 1e-10 1e-10]); */
    memset(&Qk_0[0], 0, 25U * sizeof(real_T));

    /* :  Rk = diag([5e-4 5e-4]); */
    Rk[1] = 0.0;
    Rk[2] = 0.0;
    Rk[0] = 0.0005;
    Rk[3] = 0.0005;

    /* :  Ia = x_hat(1); */
    /* :  Ib = x_hat(2); */
    /* :  Wm = x_hat(3); */
    /* :  Th_m = x_hat(4); */
    /* :  Tx = x_hat(5); */
    /* :  Th_e = Nr*Th_m; */
    Th_e = Simulacion_stepper_P.P * Simulacion_stepper_DW.x_hat[3];

    /* :  S = sin(Th_e); */
    S = sin(Th_e);

    /* :  C = cos(Th_e); */
    Th_e = cos(Th_e);

    /* :  Va = U(1); */
    /* :  Vb = U(2); */
    /* :  dx = zeros(5,1); */
    /* :  dx(1) = (Va - R*Ia + Kt*Wm*S)/L; */
    x_pred_0[0] = ((Simulacion_stepper_B.TmpSignalConversionAtSFunctio_n[0] -
                    Simulacion_stepper_P.R * Simulacion_stepper_DW.x_hat[0]) +
                   Simulacion_stepper_P.Kt * Simulacion_stepper_DW.x_hat[2] * S)
      / Simulacion_stepper_P.L;

    /* :  dx(2) = (Vb - R*Ib - Kt*Wm*C)/L; */
    x_pred_0[1] = ((Simulacion_stepper_B.TmpSignalConversionAtSFunctio_n[1] -
                    Simulacion_stepper_P.R * Simulacion_stepper_DW.x_hat[1]) -
                   Simulacion_stepper_P.Kt * Simulacion_stepper_DW.x_hat[2] *
                   Th_e) / Simulacion_stepper_P.L;

    /* :  dx(3) = (Kt*(Ib*C-Ia*S)-B*Wm-Tx)/J; */
    x_pred_0[2] = (((Simulacion_stepper_DW.x_hat[1] * Th_e -
                     Simulacion_stepper_DW.x_hat[0] * S) *
                    Simulacion_stepper_P.Kt - Simulacion_stepper_P.B *
                    Simulacion_stepper_DW.x_hat[2]) -
                   Simulacion_stepper_DW.x_hat[4]) / Simulacion_stepper_P.J;

    /* :  dx(4) = Wm; */
    x_pred_0[3] = Simulacion_stepper_DW.x_hat[2];

    /* :  dx(5) = 0; */
    x_pred_0[4] = 0.0;

    /* :  x_pred = x_hat + Ts*dx; */
    for (r1 = 0; r1 < 5; r1++) {
      Qk_0[r1 + 5 * r1] = b_0[r1];
      x_pred_0[r1] = Simulacion_stepper_P.sample_time * x_pred_0[r1] +
        Simulacion_stepper_DW.x_hat[r1];
    }

    /* :  x_pred(4) = wrap_2pi(x_pred(4)); */
    /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
    RL = Simulacion_stepper_mod(x_pred_0[3] + 3.1415926535897931);

    /* 'wrap_2pi:3' if angle < 0 */
    if (RL < 0.0) {
      /* 'wrap_2pi:4' angle = angle + 2*pi; */
      RL += 6.2831853071795862;
    }

    /* 'wrap_2pi:7' angle_output = angle - pi; */
    x_pred_0[3] = RL - 3.1415926535897931;

    /* :  RL = R/L; */
    RL = Simulacion_stepper_P.R / Simulacion_stepper_P.L;

    /* :  uno = Kt*S; */
    uno = Simulacion_stepper_P.Kt * S;

    /* :  dos = Kt*C; */
    dos = Simulacion_stepper_P.Kt * Th_e;

    /* :  A = [-RL 0   (uno/L) (dos*Nr*Wm/L) 0; */
    /* :   0  -RL   -dos/L    uno*Nr*Wm/L  0; */
    /* :   -uno/J dos/J -B/J -Kt*Nr*(Ia*C+Ib*S)/J -1/J; */
    /* :   0 0 1 0 0; */
    /* :   0 0 0 0 0]; */
    /* :  A_hat = eye(5) + Ts*A; */
    memset(&P_pred_0[0], 0, 25U * sizeof(real_T));
    tmp_0[0] = Simulacion_stepper_P.sample_time * -RL;
    tmp_0[1] = Simulacion_stepper_P.sample_time * 0.0;
    tmp_0[2] = uno / Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    tmp_0[3] = dos * Simulacion_stepper_P.P * Simulacion_stepper_DW.x_hat[2] /
      Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    tmp_0[4] = Simulacion_stepper_P.sample_time * 0.0;
    tmp_1[0] = Simulacion_stepper_P.sample_time * 0.0;
    tmp_1[1] = Simulacion_stepper_P.sample_time * -RL;
    tmp_1[2] = -dos / Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    tmp_1[3] = uno * Simulacion_stepper_P.P * Simulacion_stepper_DW.x_hat[2] /
      Simulacion_stepper_P.L * Simulacion_stepper_P.sample_time;
    tmp_1[4] = Simulacion_stepper_P.sample_time * 0.0;
    tmp_2[0] = -uno / Simulacion_stepper_P.J * Simulacion_stepper_P.sample_time;
    tmp_2[1] = dos / Simulacion_stepper_P.J * Simulacion_stepper_P.sample_time;
    tmp_2[2] = -Simulacion_stepper_P.B / Simulacion_stepper_P.J *
      Simulacion_stepper_P.sample_time;
    tmp_2[3] = (Simulacion_stepper_DW.x_hat[0] * Th_e +
                Simulacion_stepper_DW.x_hat[1] * S) * (-Simulacion_stepper_P.Kt *
      Simulacion_stepper_P.P) / Simulacion_stepper_P.J *
      Simulacion_stepper_P.sample_time;
    tmp_2[4] = -1.0 / Simulacion_stepper_P.J * Simulacion_stepper_P.sample_time;
    for (r1 = 0; r1 < 5; r1++) {
      P_pred_0[r1 + 5 * r1] = 1.0;
      tmp_d = _mm_add_pd(_mm_loadu_pd(&P_pred_0[5 * r1]), _mm_set_pd(tmp_1[r1],
        tmp_0[r1]));
      _mm_storeu_pd(&A_hat_0[5 * r1], tmp_d);
      jm1 = 5 * r1 + 2;
      A_hat_0[jm1] = P_pred_0[jm1] + tmp_2[r1];
      r2 = 5 * r1 + 3;
      tmp_d = _mm_add_pd(_mm_loadu_pd(&P_pred_0[r2]), _mm_mul_pd(_mm_set1_pd
        (Simulacion_stepper_P.sample_time), _mm_set_pd(0.0, c[r1])));
      _mm_storeu_pd(&A_hat_0[r2], tmp_d);
    }

    /* :  A_hat_T = transpose(A_hat); */
    /* :  P_pred = A_hat * P * A_hat_T + Qk; */
    for (r2 = 0; r2 < 5; r2++) {
      for (r1 = 0; r1 < 5; r1++) {
        Th_e = 0.0;
        for (jm1 = 0; jm1 < 5; jm1++) {
          Th_e += A_hat_0[5 * jm1 + r2] * Simulacion_stepper_DW.P_p[5 * r1 + jm1];
        }

        A_hat_2[r2 + 5 * r1] = Th_e;
      }

      for (r1 = 0; r1 < 5; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 5; jm1++) {
          S += A_hat_2[5 * jm1 + r2] * A_hat_0[5 * jm1 + r1];
        }

        jm1 = 5 * r1 + r2;
        P_pred_0[jm1] = Qk_0[jm1] + S;
      }
    }

    /* :  Ck = [1 0 0 0 0; */
    /* :        0 1 0 0 0]; */
    /* :  Ck_T = transpose(Ck); */
    /* :  y_pred = Ck*x_pred; */
    /* :  Sk = Ck*P_pred*Ck_T + Rk; */
    for (r2 = 0; r2 < 2; r2++) {
      for (r1 = 0; r1 < 5; r1++) {
        Th_e = 0.0;
        for (jm1 = 0; jm1 < 5; jm1++) {
          Th_e += (real_T)a_3[(jm1 << 1) + r2] * P_pred_0[5 * r1 + jm1];
        }

        a[r2 + (r1 << 1)] = Th_e;
      }

      for (r1 = 0; r1 < 2; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 5; jm1++) {
          S += a[(jm1 << 1) + r2] * (real_T)b_b_3[5 * r1 + jm1];
        }

        jm1 = (r1 << 1) + r2;
        Sk[jm1] = Rk[jm1] + S;
      }
    }

    /* :  Lk = P_pred * Ck_T /Sk; */
    for (r2 = 0; r2 < 2; r2++) {
      for (r1 = 0; r1 < 5; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 5; jm1++) {
          S += P_pred_0[5 * jm1 + r1] * (real_T)b_b_3[5 * r2 + jm1];
        }

        a[r1 + 5 * r2] = S;
      }
    }

    if (fabs(Sk[1]) > fabs(Sk[0])) {
      r1 = 1;
      r2 = 0;
    } else {
      r1 = 0;
      r2 = 1;
    }

    S = Sk[r2] / Sk[r1];
    tmp_3 = Sk[r1 + 2];
    Th_e = Sk[r2 + 2] - tmp_3 * S;
    for (jm1 = 0; jm1 < 5; jm1++) {
      Lk_tmp = 5 * r1 + jm1;
      Lk_0[Lk_tmp] = a[jm1] / Sk[r1];
      idxAjj = 5 * r2 + jm1;
      Lk_0[idxAjj] = (a[jm1 + 5] - Lk_0[Lk_tmp] * tmp_3) / Th_e;
      Lk_0[Lk_tmp] -= Lk_0[idxAjj] * S;
    }

    /* :  y_meas = [X(1);X(2)]; */
    /* :  innov = y_meas - y_pred; */
    /* :  x_corr = x_pred + Lk*(innov); */
    y_pred[0] = Simulacion_stepper_B.TmpSignalConversionAtSFunctio_h[0];
    y_pred[1] = Simulacion_stepper_B.TmpSignalConversionAtSFunctio_h[1];
    for (r2 = 0; r2 < 2; r2++) {
      Th_e = 0.0;
      for (r1 = 0; r1 < 5; r1++) {
        Th_e += (real_T)a_3[(r1 << 1) + r2] * x_pred_0[r1];
      }

      tmp[r2] = y_pred[r2] - Th_e;
    }

    S = tmp[0];
    Th_e = tmp[1];
    for (r2 = 0; r2 <= 2; r2 += 2) {
      tmp_d = _mm_loadu_pd(&Lk_0[r2]);
      tmp_d = _mm_mul_pd(tmp_d, _mm_set1_pd(S));
      tmp_b = _mm_loadu_pd(&Lk_0[r2 + 5]);
      tmp_b = _mm_mul_pd(tmp_b, _mm_set1_pd(Th_e));
      tmp_d = _mm_add_pd(tmp_b, tmp_d);
      tmp_b = _mm_loadu_pd(&x_pred_0[r2]);
      tmp_d = _mm_add_pd(tmp_b, tmp_d);
      _mm_storeu_pd(&Simulacion_stepper_B.x_corr_c[r2], tmp_d);
    }

    for (r2 = 4; r2 < 5; r2++) {
      tmp_3 = Lk_0[r2] * S;
      tmp_3 += Lk_0[r2 + 5] * Th_e;
      tmp_3 += x_pred_0[r2];
      Simulacion_stepper_B.x_corr_c[r2] = tmp_3;
    }

    /* :  x_corr(4) = wrap_2pi(x_corr(4)); */
    /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
    S = Simulacion_stepper_mod(Simulacion_stepper_B.x_corr_c[3] +
      3.1415926535897931);

    /* 'wrap_2pi:3' if angle < 0 */
    if (S < 0.0) {
      /* 'wrap_2pi:4' angle = angle + 2*pi; */
      S += 6.2831853071795862;
    }

    /* 'wrap_2pi:7' angle_output = angle - pi; */
    Simulacion_stepper_B.x_corr_c[3] = S - 3.1415926535897931;

    /* :  x_hat = x_corr; */
    /* :  Lk_T = transpose(Lk); */
    /* :  P = P_pred - Lk*Sk*Lk_T; */
    for (r1 = 0; r1 < 5; r1++) {
      Simulacion_stepper_DW.x_hat[r1] = Simulacion_stepper_B.x_corr_c[r1];
      S = Lk_0[r1];
      Th_e = S * Sk[0];
      tmp_4 = Lk_0[r1 + 5];
      Th_e += tmp_4 * Sk[1];
      tmp_3 = Th_e;
      Th_e = S * Sk[2];
      Th_e += tmp_4 * Sk[3];
      for (r2 = 0; r2 < 5; r2++) {
        S = tmp_3 * Lk_0[r2];
        S += Lk_0[r2 + 5] * Th_e;
        jm1 = 5 * r2 + r1;
        Simulacion_stepper_DW.P_p[jm1] = P_pred_0[jm1] - S;
      }
    }

    /* :  P_T = transpose(P); */
    /* :  P = 0.5 * (P + P_T); */
    for (r2 = 0; r2 < 5; r2++) {
      for (r1 = 0; r1 < 5; r1++) {
        jm1 = 5 * r2 + r1;
        Qk_0[jm1] = (Simulacion_stepper_DW.P_p[5 * r1 + r2] +
                     Simulacion_stepper_DW.P_p[jm1]) * 0.5;
      }
    }

    memcpy(&Simulacion_stepper_DW.P_p[0], &Qk_0[0], 25U * sizeof(real_T));

    /* SignalConversion generated from: '<S94>/ SFunction ' incorporates:
     *  MATLAB Function: '<S23>/MATLAB Function2'
     */
    Simulacion_stepper_B.TmpSignalConversionAtSFunctionI[0] =
      Simulacion_stepper_B.alpha;
    Simulacion_stepper_B.TmpSignalConversionAtSFunctionI[1] =
      Simulacion_stepper_B.beta;

    /* SignalConversion generated from: '<S94>/ SFunction ' incorporates:
     *  MATLAB Function: '<S23>/MATLAB Function2'
     */
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_j[0] = 0.0;
    Simulacion_stepper_B.TmpSignalConversionAtSFunctio_j[1] = 0.0;

    /* MATLAB Function: '<S23>/MATLAB Function2' */
    /* :  Nx = 4; */
    /* :  Nw = 4; */
    /* :  Ny = 2; */
    /* :  Na = Nx+Nw+Ny; */
    /* :  if isempty(x_hat) */
    /* :  x_hat = zeros(Nx,1); */
    /* :  x_corr = zeros(Nx,1); */
    /* :  x_pred = zeros(Nx,1); */
    /* :  x_hat_a = zeros(Na,1); */
    /* :  X_minus = zeros(Nx,2*Na+1); */
    /* :  Y_minus = zeros(Ny,2*Na+1); */
    /* :  x_hat_a = [x_hat;zeros(Nw,1);zeros(Ny,1)]; */
    /* :  P_a = blkdiag(P,Qk,Rk); */
    memset(&S_0[0], 0, 100U * sizeof(real_T));
    for (r2 = 0; r2 < 4; r2++) {
      Lk_tmp = r2 << 2;
      S_0[10 * r2] = Simulacion_stepper_DW.P[Lk_tmp];
      r1 = (r2 + 4) * 10;
      S_0[r1 + 4] = Simulacion_stepper_DW.Qk[Lk_tmp];
      jm1 = Lk_tmp + 1;
      S_0[10 * r2 + 1] = Simulacion_stepper_DW.P[jm1];
      S_0[r1 + 5] = Simulacion_stepper_DW.Qk[jm1];
      jm1 = Lk_tmp + 2;
      S_0[10 * r2 + 2] = Simulacion_stepper_DW.P[jm1];
      S_0[r1 + 6] = Simulacion_stepper_DW.Qk[jm1];
      Lk_tmp += 3;
      S_0[10 * r2 + 3] = Simulacion_stepper_DW.P[Lk_tmp];
      S_0[r1 + 7] = Simulacion_stepper_DW.Qk[Lk_tmp];
    }

    S_0[88] = Simulacion_stepper_DW.Rk[0];
    S_0[89] = Simulacion_stepper_DW.Rk[1];
    S_0[98] = Simulacion_stepper_DW.Rk[2];
    S_0[99] = Simulacion_stepper_DW.Rk[3];

    /* :  S = chol(P_a,"lower"); */
    r2 = 0;
    r1 = 0;
    exitg1 = false;
    while ((!exitg1) && (r1 < 10)) {
      jm1 = r1;
      idxAjj = r1 * 10 + r1;
      S = 0.0;
      if (r1 >= 1) {
        b_ix = r1;
        for (iy = 0; iy < jm1; iy++) {
          Th_e = S_0[iy * 10 + b_ix];
          S += Th_e * Th_e;
        }
      }

      S = S_0[idxAjj] - S;
      if (S > 0.0) {
        S = sqrt(S);
        S_0[idxAjj] = S;
        if (r1 + 1 < 10) {
          if (r1 != 0) {
            b_ix = ((r1 - 1) * 10 + r1) + 2;
            for (b_iy = r1 + 2; b_iy <= b_ix; b_iy += 10) {
              Lk_tmp = b_iy - r1;
              Th_e = -S_0[div_nde_s32_floor(Lk_tmp - 2, 10) * 10 + jm1];
              iy = idxAjj + 1;
              e = Lk_tmp + 8;
              for (ia = b_iy; ia <= e; ia++) {
                Lk_tmp = (iy + ia) - b_iy;
                S_0[Lk_tmp] += S_0[ia - 1] * Th_e;
              }
            }
          }

          S = 1.0 / S;
          b_ix = (idxAjj - r1) + 10;
          Lk_tmp = (((((b_ix - idxAjj) - 1) / 2) << 1) + idxAjj) + 2;
          iy = Lk_tmp - 2;
          for (jm1 = idxAjj + 2; jm1 <= iy; jm1 += 2) {
            tmp_d = _mm_loadu_pd(&S_0[jm1 - 1]);
            tmp_d = _mm_mul_pd(tmp_d, _mm_set1_pd(S));
            _mm_storeu_pd(&S_0[jm1 - 1], tmp_d);
          }

          for (jm1 = Lk_tmp; jm1 <= b_ix; jm1++) {
            S_0[jm1 - 1] *= S;
          }
        }

        r1++;
      } else {
        S_0[idxAjj] = S;
        r2 = r1 + 1;
        exitg1 = true;
      }
    }

    if (r2 == 0) {
      r2 = 11;
    }

    r2--;
    for (r1 = 2; r1 <= r2; r1++) {
      idxAjj = r1;
      memset(&S_0[r1 * 10 + -10], 0, (uint32_T)(idxAjj - 1) * sizeof(real_T));
    }

    /* :  X =  x_hat_a(:,ones(1,2*Na+1)) + h*[zeros(Na,1), S, -S]; */
    for (r2 = 0; r2 < 10; r2++) {
      X[r2] = Simulacion_stepper_DW.h * 0.0;
      for (r1 = 0; r1 <= 8; r1 += 2) {
        tmp_d = _mm_loadu_pd(&S_0[10 * r2 + r1]);
        tmp_b = _mm_mul_pd(_mm_set1_pd(Simulacion_stepper_DW.h), tmp_d);
        _mm_storeu_pd(&X[r1 + 10 * (r2 + 1)], tmp_b);
        tmp_d = _mm_mul_pd(tmp_d, _mm_set1_pd(-1.0));
        tmp_d = _mm_mul_pd(_mm_set1_pd(Simulacion_stepper_DW.h), tmp_d);
        _mm_storeu_pd(&X[r1 + 10 * (r2 + 11)], tmp_d);
      }
    }

    /* :  Va = U(1); */
    /* :  Vb = U(2); */
    /* :  for k=1:2*Na+1 */
    tmp_3 = Simulacion_stepper_B.TmpSignalConversionAtSFunctionI[0];
    tmp_4 = Simulacion_stepper_B.TmpSignalConversionAtSFunctionI[1];
    for (r1 = 0; r1 < 21; r1++) {
      /* :  xk = X(1:Nx,k); */
      /* :  wk = X(Nx+1 : Nx+Nw,k); */
      /* :  vk = X(Nx+Nw+1 : end,k); */
      /* :  Ia = xk(1); */
      /* :  Ib = xk(2); */
      /* :  Wm = xk(3); */
      /* :  Th_m = xk(4); */
      /* :  Th_e = Nr*Th_m; */
      tmp_5 = X[10 * r1 + 3];
      Th_e = 0.0 * tmp_5;

      /* :  S = sin(Th_e); */
      S = sin(Th_e);

      /* :  C = cos(Th_e); */
      Th_e = cos(Th_e);

      /* :  dx1 = (Va - R*Ia + Kt*Wm*S)/L; */
      /* :  dx2 = (Vb - R*Ib - Kt*Wm*C)/L; */
      /* :  dx3 = (Kt*(Ib*C-Ia*S)-B*Wm)/J; */
      /* :  dx4 = Wm; */
      /* :  X_minus(1,k) = Ia + Ts*dx1+ wk(1); */
      tmp_7 = X[10 * r1];
      tmp_6 = X[10 * r1 + 2];
      S = (((tmp_3 - 0.0 * tmp_7) + 0.0 * tmp_6 * S > 0.0 ? (rtInf) : (tmp_3 -
             0.0 * tmp_7) + 0.0 * tmp_6 * S < 0.0 ? (rtMinusInf) : (rtNaN)) *
           0.0 + tmp_7) + X[10 * r1 + 4];
      r2 = r1 << 2;
      X_minus[r2] = S;

      /* :  X_minus(2,k) = Ib + Ts*dx2 + wk(2); */
      tmp_7 = X[10 * r1 + 1];
      Th_e = (((tmp_4 - 0.0 * tmp_7) - 0.0 * tmp_6 * Th_e > 0.0 ? (rtInf) :
               (tmp_4 - 0.0 * tmp_7) - 0.0 * tmp_6 * Th_e < 0.0 ? (rtMinusInf) :
               (rtNaN)) * 0.0 + tmp_7) + X[10 * r1 + 5];
      X_minus[r2 + 1] = Th_e;

      /* :  X_minus(3,k) = Wm + Ts*dx3 + wk(3); */
      X_minus[r2 + 2] = (rtNaN);

      /* :  X_minus(4,k) = Th_m + Ts*dx4 + wk(4); */
      tmp_5 = (0.0 * tmp_6 + tmp_5) + X[10 * r1 + 7];

      /* :  X_minus(4,k) = wrap_2pi(X_minus(4,k)); */
      /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
      RL = Simulacion_stepper_mod(tmp_5 + 3.1415926535897931);

      /* 'wrap_2pi:3' if angle < 0 */
      if (RL < 0.0) {
        /* 'wrap_2pi:4' angle = angle + 2*pi; */
        RL += 6.2831853071795862;
      }

      /* 'wrap_2pi:7' angle_output = angle - pi; */
      tmp_5 = RL - 3.1415926535897931;
      X_minus[r2 + 3] = tmp_5;

      /* :  Y_minus(1,k) = X_minus(1,k) + vk(1); */
      jm1 = r1 << 1;
      Y_minus[jm1] = X[10 * r1 + 8] + S;

      /* :  Y_minus(2,k) = X_minus(2,k) + vk(2); */
      Y_minus[jm1 + 1] = X[10 * r1 + 9] + Th_e;
    }

    /* :  x_pred = X_minus*alpham; */
    for (r2 = 0; r2 < 4; r2++) {
      tmp_3 = 0.0;
      for (r1 = 0; r1 < 21; r1++) {
        tmp_3 += X_minus[(r1 << 2) + r2] * Simulacion_stepper_DW.alpham[r1];
      }

      x_pred[r2] = tmp_3;
    }

    /* :  y_pred = Y_minus*alpham; */
    for (r2 = 0; r2 < 2; r2++) {
      Th_e = 0.0;
      for (r1 = 0; r1 < 21; r1++) {
        Th_e += Y_minus[(r1 << 1) + r2] * Simulacion_stepper_DW.alpham[r1];
      }

      y_pred[r2] = Th_e;
    }

    /* :  Xs = X_minus - x_pred(:,ones([1 2*Na+1])); */
    /* :  Ys = Y_minus - y_pred(:,ones([1 2*Na+1])); */
    tmp_3 = x_pred[0];
    S = x_pred[1];
    tmp_4 = x_pred[2];
    tmp_5 = x_pred[3];
    Th_e = y_pred[0];
    tmp_6 = y_pred[1];
    for (r2 = 0; r2 <= 18; r2 += 2) {
      r1 = r2 << 2;
      jm1 = (r2 + 1) << 2;
      tmp_d = _mm_sub_pd(_mm_set_pd(X_minus[jm1], X_minus[r1]), _mm_set1_pd
                         (tmp_3));
      _mm_storeu_pd(&tmp_c[0], tmp_d);
      Xs[r1] = tmp_c[0];
      Xs[jm1] = tmp_c[1];
      Lk_tmp = r1 + 1;
      idxAjj = jm1 + 1;
      tmp_d = _mm_sub_pd(_mm_set_pd(X_minus[idxAjj], X_minus[Lk_tmp]),
                         _mm_set1_pd(S));
      _mm_storeu_pd(&tmp_c[0], tmp_d);
      Xs[Lk_tmp] = tmp_c[0];
      Xs[idxAjj] = tmp_c[1];
      Lk_tmp = r1 + 2;
      idxAjj = jm1 + 2;
      tmp_d = _mm_sub_pd(_mm_set_pd(X_minus[idxAjj], X_minus[Lk_tmp]),
                         _mm_set1_pd(tmp_4));
      _mm_storeu_pd(&tmp_c[0], tmp_d);
      Xs[Lk_tmp] = tmp_c[0];
      Xs[idxAjj] = tmp_c[1];
      r1 += 3;
      jm1 += 3;
      tmp_d = _mm_sub_pd(_mm_set_pd(X_minus[jm1], X_minus[r1]), _mm_set1_pd
                         (tmp_5));
      _mm_storeu_pd(&tmp_c[0], tmp_d);
      Xs[r1] = tmp_c[0];
      Xs[jm1] = tmp_c[1];
      r1 = r2 << 1;
      jm1 = (r2 + 1) << 1;
      tmp_d = _mm_sub_pd(_mm_set_pd(Y_minus[jm1], Y_minus[r1]), _mm_set1_pd(Th_e));
      _mm_storeu_pd(&tmp_c[0], tmp_d);
      Ys[r1] = tmp_c[0];
      Ys[jm1] = tmp_c[1];
      r1++;
      jm1++;
      tmp_d = _mm_sub_pd(_mm_set_pd(Y_minus[jm1], Y_minus[r1]), _mm_set1_pd
                         (tmp_6));
      _mm_storeu_pd(&tmp_c[0], tmp_d);
      Ys[r1] = tmp_c[0];
      Ys[jm1] = tmp_c[1];
    }

    for (r2 = 20; r2 < 21; r2++) {
      r1 = r2 << 2;
      tmp_d = _mm_sub_pd(_mm_loadu_pd(&X_minus[r1]), _mm_set_pd(S, tmp_3));
      _mm_storeu_pd(&Xs[r1], tmp_d);
      r1 += 2;
      tmp_d = _mm_sub_pd(_mm_loadu_pd(&X_minus[r1]), _mm_set_pd(tmp_5, tmp_4));
      _mm_storeu_pd(&Xs[r1], tmp_d);
      r1 = r2 << 1;
      tmp_d = _mm_sub_pd(_mm_loadu_pd(&Y_minus[r1]), _mm_set_pd(tmp_6, Th_e));
      _mm_storeu_pd(&Ys[r1], tmp_d);
    }

    /* :  Xs_T = transpose(Xs); */
    /* :  Ys_T = transpose(Ys); */
    for (r2 = 0; r2 < 21; r2++) {
      jm1 = r2 << 1;
      Y_minus[r2] = Ys[jm1];
      Y_minus[r2 + 21] = Ys[jm1 + 1];
    }

    /* :  P_pred = Xs*diag(alphac)*Xs_T; */
    memset(&b[0], 0, 441U * sizeof(real_T));
    for (r1 = 0; r1 < 21; r1++) {
      b[r1 + 21 * r1] = Simulacion_stepper_DW.alphac[r1];
    }

    /* :  Pyy = Ys* diag(alphac)* Ys_T; */
    memset(&b_b[0], 0, 441U * sizeof(real_T));
    for (r1 = 0; r1 < 21; r1++) {
      b_b[r1 + 21 * r1] = Simulacion_stepper_DW.alphac[r1];
    }

    for (r2 = 0; r2 < 2; r2++) {
      for (r1 = 0; r1 < 21; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 21; jm1++) {
          S += Ys[(jm1 << 1) + r2] * b_b[21 * r1 + jm1];
        }

        Ys_0[r2 + (r1 << 1)] = S;
      }

      for (r1 = 0; r1 < 2; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 21; jm1++) {
          S += Ys_0[(jm1 << 1) + r2] * Y_minus[21 * r1 + jm1];
        }

        Rk[r2 + (r1 << 1)] = S;
      }
    }

    /* :  Pxy = Xs*diag(alphac)*Ys_T; */
    memset(&b_b[0], 0, 441U * sizeof(real_T));
    for (r1 = 0; r1 < 21; r1++) {
      b_b[r1 + 21 * r1] = Simulacion_stepper_DW.alphac[r1];
    }

    for (r2 = 0; r2 < 4; r2++) {
      for (r1 = 0; r1 < 21; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 21; jm1++) {
          S += Xs[(jm1 << 2) + r2] * b_b[21 * r1 + jm1];
        }

        X_minus[r2 + (r1 << 2)] = S;
      }

      for (r1 = 0; r1 < 2; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 21; jm1++) {
          S += X_minus[(jm1 << 2) + r2] * Y_minus[21 * r1 + jm1];
        }

        B[r2 + (r1 << 2)] = S;
      }
    }

    /* :  Lk = Pxy/Pyy; */
    if (fabs(Rk[1]) > fabs(Rk[0])) {
      r1 = 1;
      r2 = 0;
    } else {
      r1 = 0;
      r2 = 1;
    }

    S = Rk[r2] / Rk[r1];
    tmp_3 = Rk[r1 + 2];
    Th_e = Rk[r2 + 2] - tmp_3 * S;
    Lk_tmp = r1 << 2;
    Lk[Lk_tmp] = B[0] / Rk[r1];
    idxAjj = r2 << 2;
    Lk[idxAjj] = (B[4] - Lk[Lk_tmp] * tmp_3) / Th_e;
    Lk[Lk_tmp] -= Lk[idxAjj] * S;
    r2 = Lk_tmp + 1;
    Lk[r2] = B[1] / Rk[r1];
    jm1 = idxAjj + 1;
    Lk[jm1] = (B[5] - Lk[r2] * tmp_3) / Th_e;
    Lk[r2] -= Lk[jm1] * S;
    r2 = Lk_tmp + 2;
    Lk[r2] = B[2] / Rk[r1];
    jm1 = idxAjj + 2;
    Lk[jm1] = (B[6] - Lk[r2] * tmp_3) / Th_e;
    Lk[r2] -= Lk[jm1] * S;
    Lk_tmp += 3;
    Lk[Lk_tmp] = B[3] / Rk[r1];
    idxAjj += 3;
    Lk[idxAjj] = (B[7] - Lk[Lk_tmp] * tmp_3) / Th_e;
    Lk[Lk_tmp] -= Lk[idxAjj] * S;

    /* :  Lk_T = transpose(Lk); */
    /* :  x_corr = x_pred + Lk*(Y - y_pred); */
    tmp_d = _mm_sub_pd(_mm_set1_pd(0.0), _mm_loadu_pd(&y_pred[0]));

    /* End of Outputs for S-Function (fcgen): '<Root>/Function-Call Generator' */
    _mm_storeu_pd(&tmp_c[0], tmp_d);

    /* S-Function (fcgen): '<Root>/Function-Call Generator' incorporates:
     *  SubSystem: '<Root>/Triggered Subsystem1'
     */
    /* MATLAB Function: '<S23>/MATLAB Function2' */
    tmp_3 = tmp_c[0];
    tmp_4 = tmp_c[1];
    for (r2 = 0; r2 <= 2; r2 += 2) {
      tmp_d = _mm_loadu_pd(&Lk[r2]);
      tmp_d = _mm_mul_pd(tmp_d, _mm_set1_pd(tmp_3));
      tmp_b = _mm_loadu_pd(&Lk[r2 + 4]);
      tmp_b = _mm_mul_pd(tmp_b, _mm_set1_pd(tmp_4));
      tmp_d = _mm_add_pd(tmp_b, tmp_d);
      tmp_b = _mm_loadu_pd(&x_pred[r2]);
      tmp_d = _mm_add_pd(tmp_b, tmp_d);
      _mm_storeu_pd(&Simulacion_stepper_B.x_corr[r2], tmp_d);
    }

    /* :  x_corr(4) = wrap_2pi(x_corr(4)); */
    /* 'wrap_2pi:2' angle = mod(angle + pi, 2*pi); */
    S = Simulacion_stepper_mod(Simulacion_stepper_B.x_corr[3] +
      3.1415926535897931);

    /* 'wrap_2pi:3' if angle < 0 */
    if (S < 0.0) {
      /* 'wrap_2pi:4' angle = angle + 2*pi; */
      S += 6.2831853071795862;
    }

    /* 'wrap_2pi:7' angle_output = angle - pi; */
    Simulacion_stepper_B.x_corr[3] = S - 3.1415926535897931;

    /* :  x_hat = x_corr; */
    /* :  P = P_pred - Lk*Pyy*Lk_T; */
    for (r2 = 0; r2 < 4; r2++) {
      for (r1 = 0; r1 < 21; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 21; jm1++) {
          S += Xs[(jm1 << 2) + r2] * b[21 * r1 + jm1];
        }

        X_minus[r2 + (r1 << 2)] = S;
      }

      S = Lk[r2];
      Th_e = S * Rk[0];
      tmp_4 = Lk[r2 + 4];
      Th_e += tmp_4 * Rk[1];
      B[r2] = Th_e;
      Th_e = S * Rk[2];
      Th_e += tmp_4 * Rk[3];
      B[r2 + 4] = Th_e;
      for (r1 = 0; r1 < 4; r1++) {
        S = 0.0;
        for (jm1 = 0; jm1 < 21; jm1++) {
          Lk_tmp = jm1 << 2;
          S += X_minus[Lk_tmp + r2] * Xs[Lk_tmp + r1];
        }

        Lk_tmp = (r1 << 2) + r2;
        Qk[Lk_tmp] = S;
        S = B[r2] * Lk[r1];
        S += B[r2 + 4] * Lk[r1 + 4];
        P_pred[Lk_tmp] = S;
      }
    }

    for (r2 = 0; r2 <= 14; r2 += 2) {
      tmp_d = _mm_loadu_pd(&Qk[r2]);
      tmp_b = _mm_loadu_pd(&P_pred[r2]);
      tmp_d = _mm_sub_pd(tmp_d, tmp_b);
      _mm_storeu_pd(&Simulacion_stepper_DW.P[r2], tmp_d);
    }

    /* End of Outputs for S-Function (fcgen): '<Root>/Function-Call Generator' */
    /* Gain: '<Root>/Gain' */
    S = 2.0 / Simulacion_stepper_P.Vdc;

    /* Gain: '<Root>/Gain' */
    Simulacion_stepper_B.Gain = S * Simulacion_stepper_B.alpha;

    /* DigitalClock: '<S91>/Digital Clock' */
    Simulacion_stepper_B.DigitalClock =
      (((Simulacion_stepper_M->Timing.clockTick1+
         Simulacion_stepper_M->Timing.clockTickH1* 4294967296.0)) * 5.0E-6);

    /* Sum: '<S91>/Add1' incorporates:
     *  Constant: '<S91>/Constant3'
     */
    Simulacion_stepper_B.Add1 = Simulacion_stepper_B.DigitalClock +
      Simulacion_stepper_P.Constant3_Value;

    /* Math: '<S91>/Math Function' incorporates:
     *  Constant: '<S91>/Constant1'
     */
    Simulacion_stepper_B.MathFunction = rt_remd_snf(Simulacion_stepper_B.Add1,
      Simulacion_stepper_P.Constant1_Value);

    /* Gain: '<S91>/1\ib1' */
    Simulacion_stepper_B.uib1 = Simulacion_stepper_P.uib1_Gain *
      Simulacion_stepper_B.MathFunction;

    /* Lookup_n-D: '<S91>/1-D Lookup Table' incorporates:
     *  Gain: '<S91>/1\ib1'
     */
    Simulacion_stepper_B.uDLookupTable = look1_pbinlxmpw
      (Simulacion_stepper_B.uib1, Simulacion_stepper_P.uDLookupTable_bp01Data,
       Simulacion_stepper_P.uDLookupTable_tableData,
       &Simulacion_stepper_DW.m_bpIndex, 2U);

    /* Sum: '<S91>/Add3' incorporates:
     *  Constant: '<S91>/Constant2'
     */
    Simulacion_stepper_B.Add3 = Simulacion_stepper_B.uDLookupTable -
      Simulacion_stepper_P.Constant2_Value;

    /* RelationalOperator: '<Root>/Relational Operator' */
    Simulacion_stepper_B.RelationalOperator = (int16_T)
      (Simulacion_stepper_B.Gain >= Simulacion_stepper_B.Add3);

    /* RelationalOperator: '<Root>/Relational Operator4' */
    Simulacion_stepper_B.RelationalOperator4 = (int16_T)
      (Simulacion_stepper_B.Gain < Simulacion_stepper_B.Add3);

    /* Sum: '<Root>/Sum6' */
    Simulacion_stepper_B.Sum6 = (int16_T)
      (Simulacion_stepper_B.RelationalOperator -
       Simulacion_stepper_B.RelationalOperator4);

    /* Gain: '<Root>/Gain1' */
    S = 2.0 / Simulacion_stepper_P.Vdc;

    /* Gain: '<Root>/Gain1' */
    Simulacion_stepper_B.Gain1 = S * Simulacion_stepper_B.beta;

    /* RelationalOperator: '<Root>/Relational Operator1' */
    Simulacion_stepper_B.RelationalOperator1 = (int16_T)
      (Simulacion_stepper_B.Gain1 >= Simulacion_stepper_B.Add3);

    /* RelationalOperator: '<Root>/Relational Operator5' */
    Simulacion_stepper_B.RelationalOperator5 = (int16_T)
      (Simulacion_stepper_B.Gain1 < Simulacion_stepper_B.Add3);

    /* Sum: '<Root>/Sum7' */
    Simulacion_stepper_B.Sum7 = (int16_T)
      (Simulacion_stepper_B.RelationalOperator1 -
       Simulacion_stepper_B.RelationalOperator5);

    /* Fcn: '<S19>/alpha' */
    Simulacion_stepper_B.alpha_p = cos(0.0 - Simulacion_stepper_B.Gain2) *
      Simulacion_stepper_B.CCCodeBlock_o3 + sin(0.0 - Simulacion_stepper_B.Gain2)
      * Simulacion_stepper_B.CCCodeBlock_o4;

    /* Fcn: '<S19>/beta' */
    Simulacion_stepper_B.beta_n = -sin(0.0 - Simulacion_stepper_B.Gain2) *
      Simulacion_stepper_B.CCCodeBlock_o3 + cos(0.0 - Simulacion_stepper_B.Gain2)
      * Simulacion_stepper_B.CCCodeBlock_o4;

    /* Gain: '<Root>/Gain3' */
    Simulacion_stepper_B.Gain3 = Simulacion_stepper_P.P *
      Simulacion_stepper_B.DiscreteTimeIntegrator1;

    /* Product: '<Root>/Product' */
    Simulacion_stepper_B.Product = Simulacion_stepper_B.CCCodeBlock_o4 *
      Simulacion_stepper_B.CCCodeBlock_o4;

    /* Product: '<Root>/Product1' */
    Simulacion_stepper_B.Product1 = Simulacion_stepper_B.CCCodeBlock_o3 *
      Simulacion_stepper_B.CCCodeBlock_o3;

    /* Sum: '<Root>/Sum13' */
    Simulacion_stepper_B.Sum13 = Simulacion_stepper_B.Product +
      Simulacion_stepper_B.Product1;

    /* Sqrt: '<Root>/Sqrt' */
    Simulacion_stepper_B.Sqrt = sqrt(Simulacion_stepper_B.Sum13);

    /* Product: '<Root>/Product2' */
    Simulacion_stepper_B.Product2 = Simulacion_stepper_B.CCCodeBlock_o2 *
      Simulacion_stepper_B.CCCodeBlock_o2;

    /* Product: '<Root>/Product3' */
    Simulacion_stepper_B.Product3 = Simulacion_stepper_B.CCCodeBlock_o1 *
      Simulacion_stepper_B.CCCodeBlock_o1;

    /* Sum: '<Root>/Sum14' */
    Simulacion_stepper_B.Sum14 = Simulacion_stepper_B.Product2 +
      Simulacion_stepper_B.Product3;

    /* Sqrt: '<Root>/Sqrt1' */
    Simulacion_stepper_B.Sqrt1 = sqrt(Simulacion_stepper_B.Sum14);

    /* Sum: '<Root>/Sum2' incorporates:
     *  Sum: '<Root>/Sum4'
     */
    tmp_d = _mm_sub_pd(_mm_set1_pd(Simulacion_stepper_B.donotdeletethisgain),
                       _mm_set_pd(Simulacion_stepper_B.x_corr_c[0],
      Simulacion_stepper_B.x_corr_a[0]));
    _mm_storeu_pd(&tmp_c[0], tmp_d);

    /* Sum: '<Root>/Sum2' */
    Simulacion_stepper_B.Sum2 = tmp_c[0];

    /* Sum: '<Root>/Sum4' */
    Simulacion_stepper_B.Sum4 = tmp_c[1];

    /* Sum: '<Root>/Sum8' */
    Simulacion_stepper_B.Sum8 = Simulacion_stepper_B.donotdeletethisgain;

    /* Gain: '<Root>/Gain7' incorporates:
     *  Gain: '<Root>/Gain6'
     */
    tmp_d = _mm_mul_pd(_mm_set1_pd(Simulacion_stepper_P.P), _mm_set_pd
                       (Simulacion_stepper_B.x_corr_c[3],
                        Simulacion_stepper_B.x_corr_a[3]));
    _mm_storeu_pd(&tmp_c[0], tmp_d);

    /* Gain: '<Root>/Gain7' */
    Simulacion_stepper_B.Gain7 = tmp_c[0];

    /* Gain: '<Root>/Gain6' */
    Simulacion_stepper_B.Gain6 = tmp_c[1];

    /* Gain: '<Root>/Gain8' */
    Simulacion_stepper_B.Gain8 = Simulacion_stepper_P.P * 0.0;

    /* Sum: '<Root>/Sum1' */
    Simulacion_stepper_B.Sum1 = Simulacion_stepper_B.Gain2 -
      Simulacion_stepper_B.Gain7;

    /* Sum: '<Root>/Sum5' */
    Simulacion_stepper_B.Sum5 = Simulacion_stepper_B.Gain2 -
      Simulacion_stepper_B.Gain6;

    /* Sum: '<Root>/Sum9' */
    Simulacion_stepper_B.Sum9 = Simulacion_stepper_B.Gain2 -
      Simulacion_stepper_B.Gain8;

    /* MATLAB Function: '<Root>/MATLAB Function3' */
    Simulacion_stepp_MATLABFunction(Simulacion_stepper_B.DiscreteTimeIntegrator,
      &Simulacion_stepper_B.sf_MATLABFunction3);

    /* MATLAB Function: '<Root>/MATLAB Function2' */
    Simulacion_step_MATLABFunction1(Simulacion_stepper_B.x_corr_a[3],
      &Simulacion_stepper_B.sf_MATLABFunction2);

    /* MATLAB Function: '<Root>/MATLAB Function1' */
    Simulacion_step_MATLABFunction1(Simulacion_stepper_B.x_corr_c[3],
      &Simulacion_stepper_B.sf_MATLABFunction1);

    /* MATLAB Function: '<Root>/MATLAB Function' */
    Simulacion_stepp_MATLABFunction(0.0, &Simulacion_stepper_B.sf_MATLABFunction);

    /* Gain: '<S88>/p' */
    Simulacion_stepper_B.p = Simulacion_stepper_P.p_Gain *
      Simulacion_stepper_B.DiscreteTimeIntegrator;

    /* Sum: '<S88>/Sum1' incorporates:
     *  Constant: '<S87>/Constant'
     */
    Simulacion_stepper_B.Sum1_k[0] = Simulacion_stepper_B.p +
      Simulacion_stepper_P.Constant_Value_f[0];

    /* Math: '<S88>/Math Function' incorporates:
     *  Constant: '<S88>/Constant1'
     */
    Simulacion_stepper_B.MathFunction_h[0] = rt_modd_snf
      (Simulacion_stepper_B.Sum1_k[0], Simulacion_stepper_P.Constant1_Value_md);

    /* Trigonometry: '<S88>/Trigonometric Function' */
    Th_e = sin(Simulacion_stepper_B.MathFunction_h[0]);
    Simulacion_stepper_B.TrigonometricFunction[0] = Th_e;

    /* Gain: '<S88>/pxPsim' */
    Th_e *= Simulacion_stepper_P.pxPsim_Gain;
    Simulacion_stepper_B.pxPsim[0] = Th_e;

    /* Product: '<S90>/Product' */
    Th_e *= Simulacion_stepper_B.IphA[0];
    Simulacion_stepper_B.Product_c[0] = Th_e;

    /* Sum: '<S90>/Sum of Elements' */
    S = Th_e;

    /* Sum: '<S88>/Sum1' incorporates:
     *  Constant: '<S87>/Constant'
     */
    Simulacion_stepper_B.Sum1_k[1] = Simulacion_stepper_B.p +
      Simulacion_stepper_P.Constant_Value_f[1];

    /* Math: '<S88>/Math Function' incorporates:
     *  Constant: '<S88>/Constant1'
     */
    Simulacion_stepper_B.MathFunction_h[1] = rt_modd_snf
      (Simulacion_stepper_B.Sum1_k[1], Simulacion_stepper_P.Constant1_Value_md);

    /* Trigonometry: '<S88>/Trigonometric Function' */
    Th_e = sin(Simulacion_stepper_B.MathFunction_h[1]);
    Simulacion_stepper_B.TrigonometricFunction[1] = Th_e;

    /* Gain: '<S88>/pxPsim' */
    Th_e *= Simulacion_stepper_P.pxPsim_Gain;
    Simulacion_stepper_B.pxPsim[1] = Th_e;

    /* Product: '<S90>/Product' */
    Th_e *= Simulacion_stepper_B.IphA[1];
    Simulacion_stepper_B.Product_c[1] = Th_e;

    /* Sum: '<S90>/Sum of Elements' */
    S += Th_e;

    /* Sum: '<S90>/Sum of Elements' */
    Simulacion_stepper_B.SumofElements = S;

    /* Gain: '<S88>/2' */
    Simulacion_stepper_B.u = Simulacion_stepper_P.u_Gain *
      Simulacion_stepper_B.p;

    /* Math: '<S88>/Math Function1' incorporates:
     *  Constant: '<S88>/Constant1'
     */
    Simulacion_stepper_B.MathFunction1 = rt_modd_snf(Simulacion_stepper_B.u,
      Simulacion_stepper_P.Constant1_Value_md);

    /* Trigonometry: '<S88>/Trigonometric Function1' */
    Simulacion_stepper_B.TrigonometricFunction1 = sin
      (Simulacion_stepper_B.MathFunction1);

    /* Gain: '<S90>/Tdm' */
    Simulacion_stepper_B.Tdm = Simulacion_stepper_P.Tdm_Gain *
      Simulacion_stepper_B.TrigonometricFunction1;

    /* Sum: '<S90>/Sum2' */
    Simulacion_stepper_B.TeNm = Simulacion_stepper_B.SumofElements +
      Simulacion_stepper_B.Tdm;

    /* DataTypeConversion: '<S36>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion =
      Simulacion_stepper_B.RelationalOperator;

    /* DataTypeConversion: '<S43>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_b =
      Simulacion_stepper_B.RelationalOperator4;

    /* Logic: '<Root>/NOT' */
    Simulacion_stepper_B.NOT = (Simulacion_stepper_B.RelationalOperator == 0);

    /* DataTypeConversion: '<S50>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_n = Simulacion_stepper_B.NOT;

    /* Logic: '<Root>/NOT4' */
    Simulacion_stepper_B.NOT4 = (Simulacion_stepper_B.RelationalOperator4 == 0);

    /* DataTypeConversion: '<S57>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_p = Simulacion_stepper_B.NOT4;

    /* DataTypeConversion: '<S64>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_k =
      Simulacion_stepper_B.RelationalOperator1;

    /* DataTypeConversion: '<S71>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_ky =
      Simulacion_stepper_B.RelationalOperator5;

    /* Logic: '<Root>/NOT1' */
    Simulacion_stepper_B.NOT1 = (Simulacion_stepper_B.RelationalOperator1 == 0);

    /* DataTypeConversion: '<S78>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_d = Simulacion_stepper_B.NOT1;

    /* Logic: '<Root>/NOT5' */
    Simulacion_stepper_B.NOT5 = (Simulacion_stepper_B.RelationalOperator5 == 0);

    /* DataTypeConversion: '<S85>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_i = Simulacion_stepper_B.NOT5;

    /* Gain: '<S89>/B' */
    Simulacion_stepper_B.B = Simulacion_stepper_P.B_Gain *
      Simulacion_stepper_B.DiscreteTimeIntegrator1;

    /* Sum: '<S89>/Sum1' incorporates:
     *  Constant: '<Root>/Constant'
     */
    Simulacion_stepper_B.Sum1_l = (Simulacion_stepper_B.TeNm -
      Simulacion_stepper_P.Constant_Value) - Simulacion_stepper_B.B;

    /* Gain: '<S89>/1//J' */
    Simulacion_stepper_B.uJ = Simulacion_stepper_P.uJ_Gain *
      Simulacion_stepper_B.Sum1_l;

    /* Product: '<S88>/Product1' */
    S = Simulacion_stepper_B.DiscreteTimeIntegrator1 *
      Simulacion_stepper_B.pxPsim[0];
    Simulacion_stepper_B.Product1_m[0] = S;

    /* Gain: '<S90>/R' */
    Th_e = Simulacion_stepper_P.R_Gain * Simulacion_stepper_B.IphA[0];
    Simulacion_stepper_B.R[0] = Th_e;

    /* Sum: '<S90>/Sum' */
    S = (Simulacion_stepper_B.StateSpace_o1[16] - S) - Th_e;
    Simulacion_stepper_B.Sum[0] = S;

    /* Gain: '<S90>/1//L' */
    Simulacion_stepper_B.IphA_e[0] = Simulacion_stepper_P.uL_Gain * S;

    /* Product: '<S88>/Product1' */
    S = Simulacion_stepper_B.DiscreteTimeIntegrator1 *
      Simulacion_stepper_B.pxPsim[1];
    Simulacion_stepper_B.Product1_m[1] = S;

    /* Gain: '<S90>/R' */
    Th_e = Simulacion_stepper_P.R_Gain * Simulacion_stepper_B.IphA[1];
    Simulacion_stepper_B.R[1] = Th_e;

    /* Sum: '<S90>/Sum' */
    S = (Simulacion_stepper_B.StateSpace_o1[17] - S) - Th_e;
    Simulacion_stepper_B.Sum[1] = S;

    /* Gain: '<S90>/1//L' */
    Simulacion_stepper_B.IphA_e[1] = Simulacion_stepper_P.uL_Gain * S;

    /* Gain: '<S25>/do not delete this gain' incorporates:
     *  Gain: '<S24>/do not delete this gain'
     */
    tmp_d = _mm_mul_pd(_mm_set_pd
                       (Simulacion_stepper_P.donotdeletethisgain_Gain_i,
                        Simulacion_stepper_P.donotdeletethisgain_Gain_j),
                       _mm_loadu_pd(&Simulacion_stepper_B.StateSpace_o1[18]));
    _mm_storeu_pd(&tmp_c[0], tmp_d);

    /* Gain: '<S24>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain_j = tmp_c[0];

    /* Gain: '<S25>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain_n = tmp_c[1];

    /* Update for DiscreteIntegrator: '<S90>/Discrete-Time Integrator' */
    tmp_d = _mm_add_pd(_mm_mul_pd(_mm_set1_pd
      (Simulacion_stepper_P.DiscreteTimeIntegrator_gainval), _mm_loadu_pd
      (&Simulacion_stepper_B.IphA_e[0])), _mm_loadu_pd
                       (&Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0]));

    /* Update for DiscreteIntegrator: '<S90>/Discrete-Time Integrator' */
    _mm_storeu_pd(&Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0], tmp_d);

    /* Update for S-Function (sfun_spssw_discc): '<S97>/State-Space' incorporates:
     *  Constant: '<S29>/DC'
     *  Constant: '<S30>/DC'
     *  Constant: '<S99>/SwitchCurrents'
     */
    {
      int_T *gState = (int_T*)Simulacion_stepper_DW.StateSpace_PWORK.G_STATE;

      /* Store switch gates values for next step */
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion;
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion_b;
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion_n;
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion_p;
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion_k;
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion_ky;
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion_d;
      *(gState++) = (int_T) Simulacion_stepper_B.DataTypeConversion_i;
      *(gState++) = (int_T) 0.0;
      *(gState++) = (int_T) 0.0;
      *(gState++) = (int_T) 0.0;
      *(gState++) = (int_T) 0.0;
      *(gState++) = (int_T) 0.0;
      *(gState++) = (int_T) 0.0;
      *(gState++) = (int_T) 0.0;
      *(gState++) = (int_T) 0.0;
    }

    /* Update for DiscreteIntegrator: '<S89>/Discrete-Time Integrator1' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator1_DSTATE +=
      Simulacion_stepper_P.DiscreteTimeIntegrator1_gainval *
      Simulacion_stepper_B.uJ;

    /* Update for DiscreteIntegrator: '<S89>/Discrete-Time Integrator' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE_o +=
      Simulacion_stepper_P.DiscreteTimeIntegrator_gainva_k *
      Simulacion_stepper_B.DiscreteTimeIntegrator1;
  }

  /* Matfile logging */
  rt_UpdateTXYLogVars(Simulacion_stepper_M->rtwLogInfo,
                      (&Simulacion_stepper_M->Timing.taskTime0));

  /* signal main to stop simulation */
  {                                    /* Sample time: [1.0E-6s, 0.0s] */
    if ((rtmGetTFinal(Simulacion_stepper_M)!=-1) &&
        !((rtmGetTFinal(Simulacion_stepper_M)-
           Simulacion_stepper_M->Timing.taskTime0) >
          Simulacion_stepper_M->Timing.taskTime0 * (DBL_EPSILON))) {
      rtmSetErrorStatus(Simulacion_stepper_M, "Simulation finished");
    }
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++Simulacion_stepper_M->Timing.clockTick0)) {
    ++Simulacion_stepper_M->Timing.clockTickH0;
  }

  Simulacion_stepper_M->Timing.taskTime0 =
    Simulacion_stepper_M->Timing.clockTick0 *
    Simulacion_stepper_M->Timing.stepSize0 +
    Simulacion_stepper_M->Timing.clockTickH0 *
    Simulacion_stepper_M->Timing.stepSize0 * 4294967296.0;
  if (Simulacion_stepper_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update absolute timer for sample time: [5.0E-6s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 5.0E-6, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    Simulacion_stepper_M->Timing.clockTick1++;
    if (!Simulacion_stepper_M->Timing.clockTick1) {
      Simulacion_stepper_M->Timing.clockTickH1++;
    }
  }

  rate_scheduler();
}

/* Model initialize function */
void Simulacion_stepper_initialize(void)
{
  /* Registration code */

  /* initialize real-time model */
  (void) memset((void *)Simulacion_stepper_M, 0,
                sizeof(RT_MODEL_Simulacion_stepper_T));
  rtmSetTFinal(Simulacion_stepper_M, 10.0);
  Simulacion_stepper_M->Timing.stepSize0 = 1.0E-6;

  /* Setup for data logging */
  {
    static RTWLogInfo rt_DataLoggingInfo;
    rt_DataLoggingInfo.loggingInterval = (NULL);
    Simulacion_stepper_M->rtwLogInfo = &rt_DataLoggingInfo;
  }

  /* Setup for data logging */
  {
    rtliSetLogXSignalInfo(Simulacion_stepper_M->rtwLogInfo, (NULL));
    rtliSetLogXSignalPtrs(Simulacion_stepper_M->rtwLogInfo, (NULL));
    rtliSetLogT(Simulacion_stepper_M->rtwLogInfo, "tout");
    rtliSetLogX(Simulacion_stepper_M->rtwLogInfo, "");
    rtliSetLogXFinal(Simulacion_stepper_M->rtwLogInfo, "");
    rtliSetLogVarNameModifier(Simulacion_stepper_M->rtwLogInfo, "rt_");
    rtliSetLogFormat(Simulacion_stepper_M->rtwLogInfo, 4);
    rtliSetLogMaxRows(Simulacion_stepper_M->rtwLogInfo, 0);
    rtliSetLogDecimation(Simulacion_stepper_M->rtwLogInfo, 1);
    rtliSetLogY(Simulacion_stepper_M->rtwLogInfo, "");
    rtliSetLogYSignalInfo(Simulacion_stepper_M->rtwLogInfo, (NULL));
    rtliSetLogYSignalPtrs(Simulacion_stepper_M->rtwLogInfo, (NULL));
  }

  /* block I/O */
  (void) memset(((void *) &Simulacion_stepper_B), 0,
                sizeof(B_Simulacion_stepper_T));

  /* states (dwork) */
  (void) memset((void *)&Simulacion_stepper_DW, 0,
                sizeof(DW_Simulacion_stepper_T));

  /* Matfile logging */
  rt_StartDataLoggingWithStartTime(Simulacion_stepper_M->rtwLogInfo, 0.0,
    rtmGetTFinal(Simulacion_stepper_M), Simulacion_stepper_M->Timing.stepSize0,
    (&rtmGetErrorStatus(Simulacion_stepper_M)));

  /* Start for S-Function (sfun_spssw_discc): '<S97>/State-Space' incorporates:
   *  Constant: '<S29>/DC'
   *  Constant: '<S30>/DC'
   *  Constant: '<S99>/SwitchCurrents'
   */

  /* S-Function block: <S97>/State-Space */
  {
    Simulacion_stepper_DW.StateSpace_PWORK.DS = (real_T*)calloc(22 * 20, sizeof
      (real_T));
    Simulacion_stepper_DW.StateSpace_PWORK.DX_COL = (real_T*)calloc(22, sizeof
      (real_T));
    Simulacion_stepper_DW.StateSpace_PWORK.TMP2 = (real_T*)calloc(20, sizeof
      (real_T));
    Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS = (int_T*)calloc(16,
      sizeof(int_T));
    Simulacion_stepper_DW.StateSpace_PWORK.SW_CHG = (int_T*)calloc(16, sizeof
      (int_T));
    Simulacion_stepper_DW.StateSpace_PWORK.G_STATE = (int_T*)calloc(16, sizeof
      (int_T));
    Simulacion_stepper_DW.StateSpace_PWORK.Y_SWITCH = (real_T*)calloc(16, sizeof
      (real_T));
    Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_TYPES = (int_T*)calloc(16,
      sizeof(int_T));
    Simulacion_stepper_DW.StateSpace_PWORK.IDX_OUT_SW = (int_T*)calloc(16,
      sizeof(int_T));
    Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS_INIT = (int_T*)calloc
      (16, sizeof(int_T));
    Simulacion_stepper_DW.StateSpace_PWORK.USWLAST = (real_T*)calloc(16, sizeof
      (real_T));
  }

  {
    real_T b;
    int32_T i;
    static const real_T b_0[5] = { 0.0001, 0.0001, 1.0E-6, 1.0E-10, 1.0E-10 };

    /* InitializeConditions for DiscreteIntegrator: '<S90>/Discrete-Time Integrator' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0] =
      Simulacion_stepper_P.DiscreteTimeIntegrator_IC;
    Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[1] =
      Simulacion_stepper_P.DiscreteTimeIntegrator_IC;

    /* InitializeConditions for S-Function (sfun_spssw_discc): '<S97>/State-Space' incorporates:
     *  Constant: '<S29>/DC'
     *  Constant: '<S30>/DC'
     *  Constant: '<S99>/SwitchCurrents'
     */
    {
      int32_T i, j;
      real_T *Ds = (real_T*)Simulacion_stepper_DW.StateSpace_PWORK.DS;

      /* Copy and transpose D */
      for (i=0; i<22; i++) {
        for (j=0; j<20; j++)
          Ds[i*20 + j] = (Simulacion_stepper_P.StateSpace_DS_param[i + j*22]);
      }

      {
        /* Switches work vectors */
        int_T *switch_status = (int_T*)
          Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS;
        int_T *gState = (int_T*)Simulacion_stepper_DW.StateSpace_PWORK.G_STATE;
        real_T *yswitch = (real_T*)
          Simulacion_stepper_DW.StateSpace_PWORK.Y_SWITCH;
        int_T *switchTypes = (int_T*)
          Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_TYPES;
        int_T *idxOutSw = (int_T*)
          Simulacion_stepper_DW.StateSpace_PWORK.IDX_OUT_SW;
        int_T *switch_status_init = (int_T*)
          Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;

        /* Initialize work vectors */
        switch_status[0] = 0;
        switch_status_init[0] = 0;
        gState[0] = (int_T) 0.0;
        yswitch[0] = 1/0.1;
        switchTypes[0] = (int_T)1.0;
        idxOutSw[0] = ((int_T)0.0) - 1;
        switch_status[1] = 0;
        switch_status_init[1] = 0;
        gState[1] = (int_T) 0.0;
        yswitch[1] = 1/0.1;
        switchTypes[1] = (int_T)1.0;
        idxOutSw[1] = ((int_T)0.0) - 1;
        switch_status[2] = 0;
        switch_status_init[2] = 0;
        gState[2] = (int_T) 0.0;
        yswitch[2] = 1/0.1;
        switchTypes[2] = (int_T)1.0;
        idxOutSw[2] = ((int_T)0.0) - 1;
        switch_status[3] = 0;
        switch_status_init[3] = 0;
        gState[3] = (int_T) 0.0;
        yswitch[3] = 1/0.1;
        switchTypes[3] = (int_T)1.0;
        idxOutSw[3] = ((int_T)0.0) - 1;
        switch_status[4] = 0;
        switch_status_init[4] = 0;
        gState[4] = (int_T) 0.0;
        yswitch[4] = 1/0.1;
        switchTypes[4] = (int_T)1.0;
        idxOutSw[4] = ((int_T)0.0) - 1;
        switch_status[5] = 0;
        switch_status_init[5] = 0;
        gState[5] = (int_T) 0.0;
        yswitch[5] = 1/0.1;
        switchTypes[5] = (int_T)1.0;
        idxOutSw[5] = ((int_T)0.0) - 1;
        switch_status[6] = 0;
        switch_status_init[6] = 0;
        gState[6] = (int_T) 0.0;
        yswitch[6] = 1/0.1;
        switchTypes[6] = (int_T)1.0;
        idxOutSw[6] = ((int_T)0.0) - 1;
        switch_status[7] = 0;
        switch_status_init[7] = 0;
        gState[7] = (int_T) 0.0;
        yswitch[7] = 1/0.1;
        switchTypes[7] = (int_T)1.0;
        idxOutSw[7] = ((int_T)0.0) - 1;
        switch_status[8] = 0;
        switch_status_init[8] = 0;
        gState[8] = (int_T) 0.0;
        yswitch[8] = 1/0.01;
        switchTypes[8] = (int_T)3.0;
        idxOutSw[8] = ((int_T)0.0) - 1;
        switch_status[9] = 0;
        switch_status_init[9] = 0;
        gState[9] = (int_T) 0.0;
        yswitch[9] = 1/0.01;
        switchTypes[9] = (int_T)3.0;
        idxOutSw[9] = ((int_T)0.0) - 1;
        switch_status[10] = 0;
        switch_status_init[10] = 0;
        gState[10] = (int_T) 0.0;
        yswitch[10] = 1/0.01;
        switchTypes[10] = (int_T)3.0;
        idxOutSw[10] = ((int_T)0.0) - 1;
        switch_status[11] = 0;
        switch_status_init[11] = 0;
        gState[11] = (int_T) 0.0;
        yswitch[11] = 1/0.01;
        switchTypes[11] = (int_T)3.0;
        idxOutSw[11] = ((int_T)0.0) - 1;
        switch_status[12] = 0;
        switch_status_init[12] = 0;
        gState[12] = (int_T) 0.0;
        yswitch[12] = 1/0.01;
        switchTypes[12] = (int_T)3.0;
        idxOutSw[12] = ((int_T)0.0) - 1;
        switch_status[13] = 0;
        switch_status_init[13] = 0;
        gState[13] = (int_T) 0.0;
        yswitch[13] = 1/0.01;
        switchTypes[13] = (int_T)3.0;
        idxOutSw[13] = ((int_T)0.0) - 1;
        switch_status[14] = 0;
        switch_status_init[14] = 0;
        gState[14] = (int_T) 0.0;
        yswitch[14] = 1/0.01;
        switchTypes[14] = (int_T)3.0;
        idxOutSw[14] = ((int_T)0.0) - 1;
        switch_status[15] = 0;
        switch_status_init[15] = 0;
        gState[15] = (int_T) 0.0;
        yswitch[15] = 1/0.01;
        switchTypes[15] = (int_T)3.0;
        idxOutSw[15] = ((int_T)0.0) - 1;
      }
    }

    /* InitializeConditions for DiscreteIntegrator: '<S89>/Discrete-Time Integrator1' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator1_DSTATE =
      Simulacion_stepper_P.DiscreteTimeIntegrator1_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S89>/Discrete-Time Integrator' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE_o =
      Simulacion_stepper_P.DiscreteTimeIntegrator_IC_l;

    /* SystemInitialize for S-Function (fcgen): '<Root>/Function-Call Generator1' incorporates:
     *  SubSystem: '<Root>/Triggered Subsystem'
     */
    /* Start for S-Function (PI_dq): '<S22>/C//C++ Code Block' incorporates:
     *  Constant: '<Root>/Constant8'
     *  Constant: '<S22>/Constant'
     *  Constant: '<S22>/Constant1'
     *  Constant: '<S22>/Constant10'
     *  Constant: '<S22>/Constant11'
     *  Constant: '<S22>/Constant2'
     *  Constant: '<S22>/Constant3'
     *  Constant: '<S22>/Constant4'
     *  Constant: '<S22>/Constant5'
     *  Constant: '<S22>/Constant6'
     *  Constant: '<S22>/Constant7'
     *  Constant: '<S22>/Constant8'
     *  Constant: '<S22>/Constant9'
     */
    b = -Simulacion_stepper_P.Kt;

    /* S-Function Block: <S22>/C//C++ Code Block */
    PI_dq_Start_wrapper();

    /* SystemInitialize for S-Function (PI_dq): '<S22>/C//C++ Code Block' incorporates:
     *  Outport: '<S22>/Out1'
     */
    Simulacion_stepper_B.CCCodeBlock_o1 = Simulacion_stepper_P.Out1_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S22>/C//C++ Code Block' incorporates:
     *  Outport: '<S22>/Out2'
     */
    Simulacion_stepper_B.CCCodeBlock_o2 = Simulacion_stepper_P.Out2_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S22>/C//C++ Code Block' incorporates:
     *  Outport: '<S22>/Out3'
     */
    Simulacion_stepper_B.CCCodeBlock_o3 = Simulacion_stepper_P.Out3_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S22>/C//C++ Code Block' incorporates:
     *  Outport: '<S22>/Out4'
     */
    Simulacion_stepper_B.CCCodeBlock_o4 = Simulacion_stepper_P.Out4_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S22>/C//C++ Code Block' incorporates:
     *  Outport: '<S22>/Out5'
     */
    Simulacion_stepper_B.CCCodeBlock_o5 = Simulacion_stepper_P.Out5_Y0;

    /* End of SystemInitialize for S-Function (fcgen): '<Root>/Function-Call Generator1' */

    /* SystemInitialize for S-Function (fcgen): '<Root>/Function-Call Generator' incorporates:
     *  SubSystem: '<Root>/Triggered Subsystem1'
     */
    /* SystemInitialize for MATLAB Function: '<S23>/MATLAB Function' */
    /* :  x_hat = zeros(4,1); */
    Simulacion_stepper_DW.x_hat_l[0] = 0.0;
    Simulacion_stepper_DW.x_hat_l[1] = 0.0;
    Simulacion_stepper_DW.x_hat_l[2] = 0.0;
    Simulacion_stepper_DW.x_hat_l[3] = 0.0;

    /* :  P = diag([1e-4 1e-4 1e-2 1e-1]); */
    memset(&Simulacion_stepper_DW.P_n[0], 0, sizeof(real_T) << 4U);
    Simulacion_stepper_DW.P_n[0] = 0.0001;
    Simulacion_stepper_DW.P_n[5] = 0.0001;
    Simulacion_stepper_DW.P_n[10] = 0.01;
    Simulacion_stepper_DW.P_n[15] = 0.1;

    /* SystemInitialize for MATLAB Function: '<S23>/MATLAB Function1' */
    /* :  x_corr = zeros(4,1); */
    /* :  x_pred = zeros(4,1); */
    /* :  x_hat = zeros(5,1); */
    for (i = 0; i < 5; i++) {
      Simulacion_stepper_DW.x_hat[i] = 0.0;
    }

    /* :  P = diag([1e-4 1e-4 1e-6 1e-10 1e-10]); */
    memset(&Simulacion_stepper_DW.P_p[0], 0, 25U * sizeof(real_T));
    for (i = 0; i < 5; i++) {
      Simulacion_stepper_DW.P_p[i + 5 * i] = b_0[i];
    }

    /* End of SystemInitialize for MATLAB Function: '<S23>/MATLAB Function1' */

    /* SystemInitialize for MATLAB Function: '<S23>/MATLAB Function2' */
    /* :  x_corr = zeros(5,1); */
    /* :  x_pred = zeros(5,1); */
    Simulacion_stepper_DW.Rk[0] = 0.0001;
    Simulacion_stepper_DW.Rk[1] = 0.0;
    Simulacion_stepper_DW.Rk[2] = 0.0;
    Simulacion_stepper_DW.Rk[3] = 0.0001;

    /* :  h =sqrt(Na); */
    Simulacion_stepper_DW.h = 3.1622776601683795;

    /* :  alpham = ones(2*Na+1,1) * (1/(2*h^2)); */
    for (i = 0; i < 21; i++) {
      Simulacion_stepper_DW.alpham[i] = 0.049999999999999989;
    }

    /* :  alpham(1) = (h^2 - Na)/(h^2); */
    Simulacion_stepper_DW.alpham[0] = 1.77635683940025E-16;

    /* :  alphac = alpham; */
    memcpy(&Simulacion_stepper_DW.alphac[0], &Simulacion_stepper_DW.alpham[0],
           21U * sizeof(real_T));

    /* :  P = diag([1e-4 1e-4 1e-6 1e-10 ]); */
    /* :  Qk = diag([1e-4 1e-4 1e-6 1e-10 ]); */
    memset(&Simulacion_stepper_DW.P[0], 0, sizeof(real_T) << 4U);
    memset(&Simulacion_stepper_DW.Qk[0], 0, sizeof(real_T) << 4U);
    Simulacion_stepper_DW.P[0] = 0.0001;
    Simulacion_stepper_DW.P[5] = 0.0001;
    Simulacion_stepper_DW.P[10] = 1.0E-6;
    Simulacion_stepper_DW.P[15] = 1.0E-10;

    /* :  Rk = [1e-4 0; 0 1e-4]; */
    Simulacion_stepper_DW.Qk[0] = 0.0001;

    /* SystemInitialize for Outport: '<S23>/ekf_4var1' */
    Simulacion_stepper_B.x_corr_a[0] = Simulacion_stepper_P.ekf_4var1_Y0;

    /* SystemInitialize for MATLAB Function: '<S23>/MATLAB Function2' */
    Simulacion_stepper_DW.Qk[5] = 0.0001;

    /* SystemInitialize for Outport: '<S23>/ekf_4var1' */
    Simulacion_stepper_B.x_corr_a[1] = Simulacion_stepper_P.ekf_4var1_Y0;

    /* SystemInitialize for MATLAB Function: '<S23>/MATLAB Function2' */
    Simulacion_stepper_DW.Qk[10] = 1.0E-6;

    /* SystemInitialize for Outport: '<S23>/ekf_4var1' */
    Simulacion_stepper_B.x_corr_a[2] = Simulacion_stepper_P.ekf_4var1_Y0;

    /* SystemInitialize for MATLAB Function: '<S23>/MATLAB Function2' */
    Simulacion_stepper_DW.Qk[15] = 1.0E-10;

    /* SystemInitialize for Outport: '<S23>/ekf_4var1' */
    Simulacion_stepper_B.x_corr_a[3] = Simulacion_stepper_P.ekf_4var1_Y0;

    /* SystemInitialize for Outport: '<S23>/ekf_5var' */
    for (i = 0; i < 5; i++) {
      Simulacion_stepper_B.x_corr_c[i] = Simulacion_stepper_P.ekf_5var_Y0;
    }

    /* End of SystemInitialize for Outport: '<S23>/ekf_5var' */
    /* End of SystemInitialize for S-Function (fcgen): '<Root>/Function-Call Generator' */
  }
}

/* Model terminate function */
void Simulacion_stepper_terminate(void)
{
  /* Terminate for S-Function (sfun_spssw_discc): '<S97>/State-Space' incorporates:
   *  Constant: '<S29>/DC'
   *  Constant: '<S30>/DC'
   *  Constant: '<S99>/SwitchCurrents'
   */

  /* S-Function block: <S97>/State-Space */
  {
    /* Free memory */
    free(Simulacion_stepper_DW.StateSpace_PWORK.DS);
    free(Simulacion_stepper_DW.StateSpace_PWORK.DX_COL);
    free(Simulacion_stepper_DW.StateSpace_PWORK.TMP2);

    /*
     * Circuit has switches*/
    free(Simulacion_stepper_DW.StateSpace_PWORK.G_STATE);
    free(Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS);
    free(Simulacion_stepper_DW.StateSpace_PWORK.SW_CHG);
    free(Simulacion_stepper_DW.StateSpace_PWORK.SWITCH_STATUS_INIT);
  }
}
