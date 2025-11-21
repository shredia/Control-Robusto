/*
 * Simulacion_stepper.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "Simulacion_stepper".
 *
 * Model version              : 1.34
 * Simulink Coder version : 24.1 (R2024a) 19-Nov-2023
 * C source code generated on : Wed Nov 19 16:44:22 2025
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Debugging
 * Validation result: Not run
 */

#include "Simulacion_stepper.h"
#include "rtwtypes.h"
#include <emmintrin.h>
#include <math.h>
#include "Simulacion_stepper_private.h"
#include <string.h>
#include "rt_nonfinite.h"
#include <float.h>

/* Block signals (default storage) */
B_Simulacion_stepper_T Simulacion_stepper_B;

/* Block states (default storage) */
DW_Simulacion_stepper_T Simulacion_stepper_DW;

/* Real-time model */
static RT_MODEL_Simulacion_stepper_T Simulacion_stepper_M_;
RT_MODEL_Simulacion_stepper_T *const Simulacion_stepper_M =
  &Simulacion_stepper_M_;
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
  __m128d tmp;
  real_T tmp_0[2];
  real_T Product1_m;
  real_T alpha_tmp;
  real_T alpha_tmp_0;
  real_T beta_tmp;
  if (Simulacion_stepper_M->Timing.TaskCounters.TID[1] == 0) {
    /* DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */
    Simulacion_stepper_B.IphA[0] =
      Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0];
    Simulacion_stepper_B.IphA[1] =
      Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[1];

    /* S-Function (sfun_spssw_discc): '<S89>/State-Space' incorporates:
     *  Constant: '<S24>/DC'
     *  Constant: '<S25>/DC'
     *  Constant: '<S91>/SwitchCurrents'
     */

    /* S-Function block: <S89>/State-Space */
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

    /* Gain: '<S82>/eliminate warning with bus selector' */
    tmp = _mm_mul_pd(_mm_set1_pd
                     (Simulacion_stepper_P.eliminatewarningwithbusselector),
                     _mm_loadu_pd(&Simulacion_stepper_B.StateSpace_o1[16]));

    /* Gain: '<S82>/eliminate warning with bus selector' */
    _mm_storeu_pd(&Simulacion_stepper_B.VphV[0], tmp);

    /* Gain: '<S2>/do not delete this gain' incorporates:
     *  Gain: '<S1>/do not delete this gain'
     */
    tmp = _mm_mul_pd(_mm_set_pd(Simulacion_stepper_P.donotdeletethisgain_Gain_l,
      Simulacion_stepper_P.donotdeletethisgain_Gain), _mm_loadu_pd
                     (&Simulacion_stepper_B.StateSpace_o1[20]));
    _mm_storeu_pd(&tmp_0[0], tmp);

    /* Gain: '<S1>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain = tmp_0[0];

    /* Gain: '<S2>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain_d = tmp_0[1];

    /* DiscreteIntegrator: '<S84>/Discrete-Time Integrator1' */
    Simulacion_stepper_B.DiscreteTimeIntegrator1 =
      Simulacion_stepper_DW.DiscreteTimeIntegrator1_DSTATE;

    /* DiscreteIntegrator: '<S84>/Discrete-Time Integrator' */
    Simulacion_stepper_B.DiscreteTimeIntegrator =
      Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE_o;

    /* Gain: '<Root>/Gain2' */
    Simulacion_stepper_B.Gain2 = Simulacion_stepper_P.P *
      Simulacion_stepper_B.DiscreteTimeIntegrator;

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
    /* S-Function (PI_dq): '<S17>/C//C++ Code Block' incorporates:
     *  Constant: '<Root>/Constant8'
     *  Constant: '<S17>/Constant'
     *  Constant: '<S17>/Constant1'
     *  Constant: '<S17>/Constant10'
     *  Constant: '<S17>/Constant11'
     *  Constant: '<S17>/Constant2'
     *  Constant: '<S17>/Constant3'
     *  Constant: '<S17>/Constant4'
     *  Constant: '<S17>/Constant5'
     *  Constant: '<S17>/Constant6'
     *  Constant: '<S17>/Constant7'
     *  Constant: '<S17>/Constant8'
     *  Constant: '<S17>/Constant9'
     */
    Product1_m = -Simulacion_stepper_P.Kt;
    PI_dq_Outputs_wrapper(&Simulacion_stepper_B.donotdeletethisgain,
                          &Simulacion_stepper_B.donotdeletethisgain_d,
                          &Simulacion_stepper_P.Constant8_Value,
                          &Simulacion_stepper_B.DiscreteTimeIntegrator1,
                          &Simulacion_stepper_P.Kp_corriente,
                          &Simulacion_stepper_P.Ki_corriente,
                          &Simulacion_stepper_P.Constant4_Value,
                          &Simulacion_stepper_P.sample_time,
                          &Simulacion_stepper_P.Vdc, &Simulacion_stepper_P.R,
                          &Simulacion_stepper_P.L, &Simulacion_stepper_P.L,
                          &Product1_m, &Simulacion_stepper_B.Gain2,
                          &Simulacion_stepper_P.Kp_w, &Simulacion_stepper_P.Ki_w,
                          &Simulacion_stepper_B.Sum3, &Simulacion_stepper_P.B,
                          &Simulacion_stepper_B.CCCodeBlock_o1,
                          &Simulacion_stepper_B.CCCodeBlock_o2,
                          &Simulacion_stepper_B.CCCodeBlock_o3,
                          &Simulacion_stepper_B.CCCodeBlock_o4,
                          &Simulacion_stepper_B.CCCodeBlock_o5);

    /* End of Outputs for S-Function (fcgen): '<Root>/Function-Call Generator1' */
    /* Fcn: '<S13>/alpha' incorporates:
     *  Fcn: '<S13>/beta'
     *  Fcn: '<S14>/alpha'
     *  Fcn: '<S14>/beta'
     */
    alpha_tmp = sin(0.0 - Simulacion_stepper_B.Gain2);
    alpha_tmp_0 = cos(0.0 - Simulacion_stepper_B.Gain2);

    /* Fcn: '<S13>/alpha' */
    Simulacion_stepper_B.alpha = alpha_tmp_0 *
      Simulacion_stepper_B.CCCodeBlock_o1 + alpha_tmp *
      Simulacion_stepper_B.CCCodeBlock_o2;

    /* Fcn: '<S13>/beta' incorporates:
     *  Fcn: '<S14>/beta'
     */
    beta_tmp = -alpha_tmp;

    /* Fcn: '<S13>/beta' */
    Simulacion_stepper_B.beta = beta_tmp * Simulacion_stepper_B.CCCodeBlock_o1 +
      alpha_tmp_0 * Simulacion_stepper_B.CCCodeBlock_o2;

    /* S-Function (ekf_6steps): '<S18>/C//C++ Code Block1' incorporates:
     *  Constant: '<S18>/Constant1'
     *  Constant: '<S18>/Constant2'
     *  Constant: '<S18>/Constant3'
     *  Constant: '<S18>/Constant4'
     *  Constant: '<S18>/Constant5'
     *  Constant: '<S18>/Constant6'
     */
    ekf_6steps_Outputs_wrapper(&Simulacion_stepper_B.alpha,
      &Simulacion_stepper_B.beta, &Simulacion_stepper_P.R,
      &Simulacion_stepper_P.L, &Simulacion_stepper_P.Kt, &Simulacion_stepper_P.J,
      &Simulacion_stepper_P.sample_time,
      &Simulacion_stepper_B.donotdeletethisgain,
      &Simulacion_stepper_B.donotdeletethisgain_d, &Simulacion_stepper_P.P,
      &Simulacion_stepper_B.CCCodeBlock1[0]);

    /* S-Function (ekf_6steps_4vars): '<S18>/C//C++ Code Block2' incorporates:
     *  Constant: '<S18>/Constant10'
     *  Constant: '<S18>/Constant11'
     *  Constant: '<S18>/Constant12'
     *  Constant: '<S18>/Constant7'
     *  Constant: '<S18>/Constant8'
     *  Constant: '<S18>/Constant9'
     */
    ekf_6steps_4vars_Outputs_wrapper(&Simulacion_stepper_B.alpha,
      &Simulacion_stepper_B.beta, &Simulacion_stepper_P.R,
      &Simulacion_stepper_P.L, &Simulacion_stepper_P.Kt, &Simulacion_stepper_P.J,
      &Simulacion_stepper_P.sample_time,
      &Simulacion_stepper_B.donotdeletethisgain,
      &Simulacion_stepper_B.donotdeletethisgain_d, &Simulacion_stepper_P.P,
      &Simulacion_stepper_B.CCCodeBlock2[0]);

    /* S-Function (SPKF): '<S18>/C//C++ Code Block3' incorporates:
     *  Constant: '<S18>/Constant13'
     *  Constant: '<S18>/Constant14'
     *  Constant: '<S18>/Constant15'
     *  Constant: '<S18>/Constant16'
     *  Constant: '<S18>/Constant17'
     *  Constant: '<S18>/Constant18'
     */
    SPKF_Outputs_wrapper(&Simulacion_stepper_B.alpha, &Simulacion_stepper_B.beta,
                         &Simulacion_stepper_P.R, &Simulacion_stepper_P.L,
                         &Simulacion_stepper_P.Kt, &Simulacion_stepper_P.J,
                         &Simulacion_stepper_P.sample_time,
                         &Simulacion_stepper_B.donotdeletethisgain,
                         &Simulacion_stepper_B.donotdeletethisgain_d,
                         &Simulacion_stepper_P.P,
                         &Simulacion_stepper_B.CCCodeBlock3_o1[0],
                         &Simulacion_stepper_B.CCCodeBlock3_o2[0]);

    /* End of Outputs for S-Function (fcgen): '<Root>/Function-Call Generator' */
    /* Gain: '<Root>/Gain' incorporates:
     *  Gain: '<Root>/Gain1'
     */
    Product1_m = 2.0 / Simulacion_stepper_P.Vdc;

    /* Gain: '<Root>/Gain' */
    Simulacion_stepper_B.Gain = Product1_m * Simulacion_stepper_B.alpha;

    /* DigitalClock: '<S86>/Digital Clock' */
    Simulacion_stepper_B.DigitalClock =
      (((Simulacion_stepper_M->Timing.clockTick1+
         Simulacion_stepper_M->Timing.clockTickH1* 4294967296.0)) * 5.0E-6);

    /* Sum: '<S86>/Add1' incorporates:
     *  Constant: '<S86>/Constant3'
     */
    Simulacion_stepper_B.Add1 = Simulacion_stepper_B.DigitalClock +
      Simulacion_stepper_P.Constant3_Value;

    /* Math: '<S86>/Math Function' incorporates:
     *  Constant: '<S86>/Constant1'
     */
    Simulacion_stepper_B.MathFunction = rt_remd_snf(Simulacion_stepper_B.Add1,
      Simulacion_stepper_P.Constant1_Value);

    /* Gain: '<S86>/1\ib1' */
    Simulacion_stepper_B.uib1 = Simulacion_stepper_P.uib1_Gain *
      Simulacion_stepper_B.MathFunction;

    /* Lookup_n-D: '<S86>/1-D Lookup Table' incorporates:
     *  Gain: '<S86>/1\ib1'
     */
    Simulacion_stepper_B.uDLookupTable = look1_pbinlxmpw
      (Simulacion_stepper_B.uib1, Simulacion_stepper_P.uDLookupTable_bp01Data,
       Simulacion_stepper_P.uDLookupTable_tableData,
       &Simulacion_stepper_DW.m_bpIndex, 2U);

    /* Sum: '<S86>/Add3' incorporates:
     *  Constant: '<S86>/Constant2'
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
    Simulacion_stepper_B.Gain1 = Product1_m * Simulacion_stepper_B.beta;

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

    /* Fcn: '<S14>/alpha' */
    Simulacion_stepper_B.alpha_p = alpha_tmp_0 *
      Simulacion_stepper_B.CCCodeBlock_o3 + alpha_tmp *
      Simulacion_stepper_B.CCCodeBlock_o4;

    /* Fcn: '<S14>/beta' */
    Simulacion_stepper_B.beta_n = beta_tmp * Simulacion_stepper_B.CCCodeBlock_o3
      + alpha_tmp_0 * Simulacion_stepper_B.CCCodeBlock_o4;

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
    tmp = _mm_sub_pd(_mm_set1_pd(Simulacion_stepper_B.donotdeletethisgain),
                     _mm_set_pd(Simulacion_stepper_B.CCCodeBlock1[0],
      Simulacion_stepper_B.CCCodeBlock2[0]));
    _mm_storeu_pd(&tmp_0[0], tmp);

    /* Sum: '<Root>/Sum2' */
    Simulacion_stepper_B.Sum2 = tmp_0[0];

    /* Sum: '<Root>/Sum4' */
    Simulacion_stepper_B.Sum4 = tmp_0[1];

    /* Sum: '<Root>/Sum8' */
    Simulacion_stepper_B.Sum8 = Simulacion_stepper_B.donotdeletethisgain -
      Simulacion_stepper_B.CCCodeBlock3_o1[0];

    /* Gain: '<Root>/Gain7' incorporates:
     *  Gain: '<Root>/Gain6'
     */
    tmp = _mm_mul_pd(_mm_set1_pd(Simulacion_stepper_P.P), _mm_set_pd
                     (Simulacion_stepper_B.CCCodeBlock1[3],
                      Simulacion_stepper_B.CCCodeBlock2[3]));
    _mm_storeu_pd(&tmp_0[0], tmp);

    /* Gain: '<Root>/Gain7' */
    Simulacion_stepper_B.Gain7 = tmp_0[0];

    /* Gain: '<Root>/Gain6' */
    Simulacion_stepper_B.Gain6 = tmp_0[1];

    /* Sum: '<Root>/Sum1' */
    Simulacion_stepper_B.Sum1 = Simulacion_stepper_B.Gain2 -
      Simulacion_stepper_B.Gain7;

    /* Sum: '<Root>/Sum5' */
    Simulacion_stepper_B.Sum5 = Simulacion_stepper_B.Gain2 -
      Simulacion_stepper_B.Gain6;

    /* Sum: '<Root>/Sum9' */
    Simulacion_stepper_B.Sum9 = Simulacion_stepper_B.Gain2 -
      Simulacion_stepper_B.CCCodeBlock3_o1[3];

    /* Math: '<Root>/Mod' incorporates:
     *  Constant: '<Root>/Constant6'
     */
    Simulacion_stepper_B.Mod = rt_modd_snf
      (Simulacion_stepper_B.DiscreteTimeIntegrator,
       Simulacion_stepper_P.Constant6_Value);

    /* Math: '<Root>/Mod2' incorporates:
     *  Constant: '<Root>/Constant1'
     */
    Simulacion_stepper_B.Mod2 = rt_modd_snf(Simulacion_stepper_B.CCCodeBlock2[3],
      Simulacion_stepper_P.Constant1_Value_j);

    /* Math: '<Root>/Mod3' incorporates:
     *  Constant: '<Root>/Constant2'
     */
    Simulacion_stepper_B.Mod3 = rt_modd_snf(Simulacion_stepper_B.CCCodeBlock1[3],
      Simulacion_stepper_P.Constant2_Value_m);

    /* Math: '<Root>/Mod1' incorporates:
     *  Constant: '<Root>/Constant3'
     */
    Simulacion_stepper_B.Mod1 = rt_modd_snf
      (Simulacion_stepper_B.CCCodeBlock3_o1[3],
       Simulacion_stepper_P.Constant3_Value_j);

    /* Gain: '<Root>/Gain9' */
    Product1_m = 1.0 / Simulacion_stepper_P.P;

    /* Gain: '<Root>/Gain9' */
    Simulacion_stepper_B.Gain9 = Product1_m * Simulacion_stepper_B.Mod1;

    /* Gain: '<S83>/p' */
    Simulacion_stepper_B.p = Simulacion_stepper_P.p_Gain *
      Simulacion_stepper_B.DiscreteTimeIntegrator;

    /* Sum: '<S83>/Sum1' incorporates:
     *  Constant: '<S82>/Constant'
     */
    Simulacion_stepper_B.Sum1_k[0] = Simulacion_stepper_B.p +
      Simulacion_stepper_P.Constant_Value_f[0];

    /* Math: '<S83>/Math Function' incorporates:
     *  Constant: '<S83>/Constant1'
     */
    Simulacion_stepper_B.MathFunction_h[0] = rt_modd_snf
      (Simulacion_stepper_B.Sum1_k[0], Simulacion_stepper_P.Constant1_Value_m);

    /* Trigonometry: '<S83>/Trigonometric Function' */
    alpha_tmp = sin(Simulacion_stepper_B.MathFunction_h[0]);
    Simulacion_stepper_B.TrigonometricFunction[0] = alpha_tmp;

    /* Gain: '<S83>/pxPsim' */
    alpha_tmp *= Simulacion_stepper_P.pxPsim_Gain;
    Simulacion_stepper_B.pxPsim[0] = alpha_tmp;

    /* Product: '<S85>/Product' */
    alpha_tmp *= Simulacion_stepper_B.IphA[0];
    Simulacion_stepper_B.Product_c[0] = alpha_tmp;

    /* Sum: '<S85>/Sum of Elements' */
    Product1_m = alpha_tmp;

    /* Sum: '<S83>/Sum1' incorporates:
     *  Constant: '<S82>/Constant'
     */
    Simulacion_stepper_B.Sum1_k[1] = Simulacion_stepper_B.p +
      Simulacion_stepper_P.Constant_Value_f[1];

    /* Math: '<S83>/Math Function' incorporates:
     *  Constant: '<S83>/Constant1'
     */
    Simulacion_stepper_B.MathFunction_h[1] = rt_modd_snf
      (Simulacion_stepper_B.Sum1_k[1], Simulacion_stepper_P.Constant1_Value_m);

    /* Trigonometry: '<S83>/Trigonometric Function' */
    alpha_tmp = sin(Simulacion_stepper_B.MathFunction_h[1]);
    Simulacion_stepper_B.TrigonometricFunction[1] = alpha_tmp;

    /* Gain: '<S83>/pxPsim' */
    alpha_tmp *= Simulacion_stepper_P.pxPsim_Gain;
    Simulacion_stepper_B.pxPsim[1] = alpha_tmp;

    /* Product: '<S85>/Product' */
    alpha_tmp *= Simulacion_stepper_B.IphA[1];
    Simulacion_stepper_B.Product_c[1] = alpha_tmp;

    /* Sum: '<S85>/Sum of Elements' */
    Product1_m += alpha_tmp;

    /* Sum: '<S85>/Sum of Elements' */
    Simulacion_stepper_B.SumofElements = Product1_m;

    /* Gain: '<S83>/2' */
    Simulacion_stepper_B.u = Simulacion_stepper_P.u_Gain *
      Simulacion_stepper_B.p;

    /* Math: '<S83>/Math Function1' incorporates:
     *  Constant: '<S83>/Constant1'
     */
    Simulacion_stepper_B.MathFunction1 = rt_modd_snf(Simulacion_stepper_B.u,
      Simulacion_stepper_P.Constant1_Value_m);

    /* Trigonometry: '<S83>/Trigonometric Function1' */
    Simulacion_stepper_B.TrigonometricFunction1 = sin
      (Simulacion_stepper_B.MathFunction1);

    /* Gain: '<S85>/Tdm' */
    Simulacion_stepper_B.Tdm = Simulacion_stepper_P.Tdm_Gain *
      Simulacion_stepper_B.TrigonometricFunction1;

    /* Sum: '<S85>/Sum2' */
    Simulacion_stepper_B.TeNm = Simulacion_stepper_B.SumofElements +
      Simulacion_stepper_B.Tdm;

    /* DataTypeConversion: '<S31>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion =
      Simulacion_stepper_B.RelationalOperator;

    /* DataTypeConversion: '<S38>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_b =
      Simulacion_stepper_B.RelationalOperator4;

    /* Logic: '<Root>/NOT' */
    Simulacion_stepper_B.NOT = (Simulacion_stepper_B.RelationalOperator == 0);

    /* DataTypeConversion: '<S45>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_n = Simulacion_stepper_B.NOT;

    /* Logic: '<Root>/NOT4' */
    Simulacion_stepper_B.NOT4 = (Simulacion_stepper_B.RelationalOperator4 == 0);

    /* DataTypeConversion: '<S52>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_p = Simulacion_stepper_B.NOT4;

    /* DataTypeConversion: '<S59>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_k =
      Simulacion_stepper_B.RelationalOperator1;

    /* DataTypeConversion: '<S66>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_ky =
      Simulacion_stepper_B.RelationalOperator5;

    /* Logic: '<Root>/NOT1' */
    Simulacion_stepper_B.NOT1 = (Simulacion_stepper_B.RelationalOperator1 == 0);

    /* DataTypeConversion: '<S73>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_d = Simulacion_stepper_B.NOT1;

    /* Logic: '<Root>/NOT5' */
    Simulacion_stepper_B.NOT5 = (Simulacion_stepper_B.RelationalOperator5 == 0);

    /* DataTypeConversion: '<S80>/Data Type Conversion' */
    Simulacion_stepper_B.DataTypeConversion_i = Simulacion_stepper_B.NOT5;

    /* Gain: '<S84>/B' */
    Simulacion_stepper_B.B = Simulacion_stepper_P.B_Gain *
      Simulacion_stepper_B.DiscreteTimeIntegrator1;

    /* Sum: '<S84>/Sum1' incorporates:
     *  Constant: '<Root>/Constant'
     */
    Simulacion_stepper_B.Sum1_l = (Simulacion_stepper_B.TeNm -
      Simulacion_stepper_P.Constant_Value) - Simulacion_stepper_B.B;

    /* Gain: '<S84>/1//J' */
    Simulacion_stepper_B.uJ = Simulacion_stepper_P.uJ_Gain *
      Simulacion_stepper_B.Sum1_l;

    /* Product: '<S83>/Product1' */
    Product1_m = Simulacion_stepper_B.DiscreteTimeIntegrator1 *
      Simulacion_stepper_B.pxPsim[0];
    Simulacion_stepper_B.Product1_m[0] = Product1_m;

    /* Gain: '<S85>/R' */
    alpha_tmp = Simulacion_stepper_P.R_Gain * Simulacion_stepper_B.IphA[0];
    Simulacion_stepper_B.R[0] = alpha_tmp;

    /* Sum: '<S85>/Sum' */
    Product1_m = (Simulacion_stepper_B.StateSpace_o1[16] - Product1_m) -
      alpha_tmp;
    Simulacion_stepper_B.Sum[0] = Product1_m;

    /* Gain: '<S85>/1//L' */
    Simulacion_stepper_B.IphA_e[0] = Simulacion_stepper_P.uL_Gain * Product1_m;

    /* Product: '<S83>/Product1' */
    Product1_m = Simulacion_stepper_B.DiscreteTimeIntegrator1 *
      Simulacion_stepper_B.pxPsim[1];
    Simulacion_stepper_B.Product1_m[1] = Product1_m;

    /* Gain: '<S85>/R' */
    alpha_tmp = Simulacion_stepper_P.R_Gain * Simulacion_stepper_B.IphA[1];
    Simulacion_stepper_B.R[1] = alpha_tmp;

    /* Sum: '<S85>/Sum' */
    Product1_m = (Simulacion_stepper_B.StateSpace_o1[17] - Product1_m) -
      alpha_tmp;
    Simulacion_stepper_B.Sum[1] = Product1_m;

    /* Gain: '<S85>/1//L' */
    Simulacion_stepper_B.IphA_e[1] = Simulacion_stepper_P.uL_Gain * Product1_m;

    /* Gain: '<S20>/do not delete this gain' incorporates:
     *  Gain: '<S19>/do not delete this gain'
     */
    tmp = _mm_mul_pd(_mm_set_pd(Simulacion_stepper_P.donotdeletethisgain_Gain_i,
      Simulacion_stepper_P.donotdeletethisgain_Gain_j), _mm_loadu_pd
                     (&Simulacion_stepper_B.StateSpace_o1[18]));
    _mm_storeu_pd(&tmp_0[0], tmp);

    /* Gain: '<S19>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain_j = tmp_0[0];

    /* Gain: '<S20>/do not delete this gain' */
    Simulacion_stepper_B.donotdeletethisgain_n = tmp_0[1];

    /* Update for DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */
    tmp = _mm_add_pd(_mm_mul_pd(_mm_set1_pd
      (Simulacion_stepper_P.DiscreteTimeIntegrator_gainval), _mm_loadu_pd
      (&Simulacion_stepper_B.IphA_e[0])), _mm_loadu_pd
                     (&Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0]));

    /* Update for DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */
    _mm_storeu_pd(&Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0], tmp);

    /* Update for S-Function (sfun_spssw_discc): '<S89>/State-Space' incorporates:
     *  Constant: '<S24>/DC'
     *  Constant: '<S25>/DC'
     *  Constant: '<S91>/SwitchCurrents'
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

    /* Update for DiscreteIntegrator: '<S84>/Discrete-Time Integrator1' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator1_DSTATE +=
      Simulacion_stepper_P.DiscreteTimeIntegrator1_gainval *
      Simulacion_stepper_B.uJ;

    /* Update for DiscreteIntegrator: '<S84>/Discrete-Time Integrator' */
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
  rtmSetTFinal(Simulacion_stepper_M, 1.0);
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

  /* Start for S-Function (sfun_spssw_discc): '<S89>/State-Space' incorporates:
   *  Constant: '<S24>/DC'
   *  Constant: '<S25>/DC'
   *  Constant: '<S91>/SwitchCurrents'
   */

  /* S-Function block: <S89>/State-Space */
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
    real_T tmp;
    int32_T i;

    /* InitializeConditions for DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[0] =
      Simulacion_stepper_P.DiscreteTimeIntegrator_IC;
    Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE[1] =
      Simulacion_stepper_P.DiscreteTimeIntegrator_IC;

    /* InitializeConditions for S-Function (sfun_spssw_discc): '<S89>/State-Space' incorporates:
     *  Constant: '<S24>/DC'
     *  Constant: '<S25>/DC'
     *  Constant: '<S91>/SwitchCurrents'
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

    /* InitializeConditions for DiscreteIntegrator: '<S84>/Discrete-Time Integrator1' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator1_DSTATE =
      Simulacion_stepper_P.DiscreteTimeIntegrator1_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S84>/Discrete-Time Integrator' */
    Simulacion_stepper_DW.DiscreteTimeIntegrator_DSTATE_o =
      Simulacion_stepper_P.DiscreteTimeIntegrator_IC_l;

    /* SystemInitialize for S-Function (fcgen): '<Root>/Function-Call Generator1' incorporates:
     *  SubSystem: '<Root>/Triggered Subsystem'
     */
    /* Start for S-Function (PI_dq): '<S17>/C//C++ Code Block' incorporates:
     *  Constant: '<Root>/Constant8'
     *  Constant: '<S17>/Constant'
     *  Constant: '<S17>/Constant1'
     *  Constant: '<S17>/Constant10'
     *  Constant: '<S17>/Constant11'
     *  Constant: '<S17>/Constant2'
     *  Constant: '<S17>/Constant3'
     *  Constant: '<S17>/Constant4'
     *  Constant: '<S17>/Constant5'
     *  Constant: '<S17>/Constant6'
     *  Constant: '<S17>/Constant7'
     *  Constant: '<S17>/Constant8'
     *  Constant: '<S17>/Constant9'
     */
    tmp = -Simulacion_stepper_P.Kt;

    /* S-Function Block: <S17>/C//C++ Code Block */
    PI_dq_Start_wrapper();

    /* SystemInitialize for S-Function (PI_dq): '<S17>/C//C++ Code Block' incorporates:
     *  Outport: '<S17>/Out1'
     */
    Simulacion_stepper_B.CCCodeBlock_o1 = Simulacion_stepper_P.Out1_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S17>/C//C++ Code Block' incorporates:
     *  Outport: '<S17>/Out2'
     */
    Simulacion_stepper_B.CCCodeBlock_o2 = Simulacion_stepper_P.Out2_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S17>/C//C++ Code Block' incorporates:
     *  Outport: '<S17>/Out3'
     */
    Simulacion_stepper_B.CCCodeBlock_o3 = Simulacion_stepper_P.Out3_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S17>/C//C++ Code Block' incorporates:
     *  Outport: '<S17>/Out4'
     */
    Simulacion_stepper_B.CCCodeBlock_o4 = Simulacion_stepper_P.Out4_Y0;

    /* SystemInitialize for S-Function (PI_dq): '<S17>/C//C++ Code Block' incorporates:
     *  Outport: '<S17>/Out5'
     */
    Simulacion_stepper_B.CCCodeBlock_o5 = Simulacion_stepper_P.Out5_Y0;

    /* End of SystemInitialize for S-Function (fcgen): '<Root>/Function-Call Generator1' */
    /* Start for S-Function (ekf_6steps): '<S18>/C//C++ Code Block1' incorporates:
     *  Constant: '<S18>/Constant1'
     *  Constant: '<S18>/Constant2'
     *  Constant: '<S18>/Constant3'
     *  Constant: '<S18>/Constant4'
     *  Constant: '<S18>/Constant5'
     *  Constant: '<S18>/Constant6'
     */

    /* S-Function Block: <S18>/C//C++ Code Block1 */
    ekf_6steps_Start_wrapper();

    /* Start for S-Function (ekf_6steps_4vars): '<S18>/C//C++ Code Block2' incorporates:
     *  Constant: '<S18>/Constant10'
     *  Constant: '<S18>/Constant11'
     *  Constant: '<S18>/Constant12'
     *  Constant: '<S18>/Constant7'
     *  Constant: '<S18>/Constant8'
     *  Constant: '<S18>/Constant9'
     */

    /* S-Function Block: <S18>/C//C++ Code Block2 */
    ekf_6steps_4vars_Start_wrapper();

    /* Start for S-Function (SPKF): '<S18>/C//C++ Code Block3' incorporates:
     *  Constant: '<S18>/Constant13'
     *  Constant: '<S18>/Constant14'
     *  Constant: '<S18>/Constant15'
     *  Constant: '<S18>/Constant16'
     *  Constant: '<S18>/Constant17'
     *  Constant: '<S18>/Constant18'
     */

    /* S-Function Block: <S18>/C//C++ Code Block3 */
    SPKF_Start_wrapper();

    /* SystemInitialize for S-Function (SPKF): '<S18>/C//C++ Code Block3' incorporates:
     *  Outport: '<S18>/X_SPKF'
     */
    Simulacion_stepper_B.CCCodeBlock3_o1[0] = Simulacion_stepper_P.X_SPKF_Y0;
    Simulacion_stepper_B.CCCodeBlock3_o1[1] = Simulacion_stepper_P.X_SPKF_Y0;
    Simulacion_stepper_B.CCCodeBlock3_o1[2] = Simulacion_stepper_P.X_SPKF_Y0;
    Simulacion_stepper_B.CCCodeBlock3_o1[3] = Simulacion_stepper_P.X_SPKF_Y0;
    for (i = 0; i < 5; i++) {
      /* SystemInitialize for S-Function (ekf_6steps): '<S18>/C//C++ Code Block1' incorporates:
       *  Outport: '<S18>/ekf_5var'
       */
      Simulacion_stepper_B.CCCodeBlock1[i] = Simulacion_stepper_P.ekf_5var_Y0;
    }

    /* SystemInitialize for S-Function (ekf_6steps_4vars): '<S18>/C//C++ Code Block2' incorporates:
     *  Outport: '<S18>/ekf_4var'
     */
    Simulacion_stepper_B.CCCodeBlock2[0] = Simulacion_stepper_P.ekf_4var_Y0;
    Simulacion_stepper_B.CCCodeBlock2[1] = Simulacion_stepper_P.ekf_4var_Y0;
    Simulacion_stepper_B.CCCodeBlock2[2] = Simulacion_stepper_P.ekf_4var_Y0;
    Simulacion_stepper_B.CCCodeBlock2[3] = Simulacion_stepper_P.ekf_4var_Y0;

    /* End of SystemInitialize for S-Function (fcgen): '<Root>/Function-Call Generator' */
  }
}

/* Model terminate function */
void Simulacion_stepper_terminate(void)
{
  /* Terminate for S-Function (sfun_spssw_discc): '<S89>/State-Space' incorporates:
   *  Constant: '<S24>/DC'
   *  Constant: '<S25>/DC'
   *  Constant: '<S91>/SwitchCurrents'
   */

  /* S-Function block: <S89>/State-Space */
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
