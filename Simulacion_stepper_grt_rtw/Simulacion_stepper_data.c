/*
 * Simulacion_stepper_data.c
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

/* Block parameters (default storage) */
P_Simulacion_stepper_T Simulacion_stepper_P = {
  /* Variable: B
   * Referenced by:
   *   '<S22>/Constant'
   *   '<S23>/Constant'
   *   '<S23>/Constant25'
   */
  0.0005,

  /* Variable: J
   * Referenced by:
   *   '<S23>/Constant22'
   *   '<S23>/Constant29'
   */
  4.7E-6,

  /* Variable: Ki_corriente
   * Referenced by: '<S22>/Constant3'
   */
  192666.66666666663,

  /* Variable: Ki_w
   * Referenced by: '<S22>/Constant11'
   */
  182936.02693602687,

  /* Variable: Kp_corriente
   * Referenced by: '<S22>/Constant2'
   */
  44.675999999999995,

  /* Variable: Kp_w
   * Referenced by: '<S22>/Constant1'
   */
  2.2823959595959589,

  /* Variable: Kt
   * Referenced by:
   *   '<S22>/Constant10'
   *   '<S23>/Constant20'
   *   '<S23>/Constant27'
   */
  0.33,

  /* Variable: L
   * Referenced by:
   *   '<S22>/Constant8'
   *   '<S22>/Constant9'
   *   '<S23>/Constant23'
   *   '<S23>/Constant30'
   */
  0.006,

  /* Variable: P
   * Referenced by:
   *   '<Root>/Gain2'
   *   '<Root>/Gain3'
   *   '<Root>/Gain6'
   *   '<Root>/Gain7'
   *   '<Root>/Gain8'
   *   '<S23>/Constant24'
   *   '<S23>/Constant31'
   */
  50.0,

  /* Variable: R
   * Referenced by:
   *   '<S22>/Constant7'
   *   '<S23>/Constant21'
   *   '<S23>/Constant28'
   */
  3.4,

  /* Variable: Tdm
   * Referenced by: '<Root>/Gain4'
   */
  0.0,

  /* Variable: Vdc
   * Referenced by:
   *   '<Root>/Gain'
   *   '<Root>/Gain1'
   *   '<S22>/Constant6'
   *   '<S29>/DC'
   *   '<S30>/DC'
   */
  24.0,

  /* Variable: sample_time
   * Referenced by:
   *   '<S22>/Constant5'
   *   '<S23>/Constant19'
   *   '<S23>/Constant26'
   */
  5.0E-6,

  /* Computed Parameter: ekf_4var1_Y0
   * Referenced by: '<S23>/ekf_4var1'
   */
  0.0,

  /* Computed Parameter: ekf_5var_Y0
   * Referenced by: '<S23>/ekf_5var'
   */
  0.0,

  /* Computed Parameter: Out1_Y0
   * Referenced by: '<S22>/Out1'
   */
  0.0,

  /* Computed Parameter: Out2_Y0
   * Referenced by: '<S22>/Out2'
   */
  0.0,

  /* Computed Parameter: Out3_Y0
   * Referenced by: '<S22>/Out3'
   */
  0.0,

  /* Computed Parameter: Out4_Y0
   * Referenced by: '<S22>/Out4'
   */
  0.0,

  /* Computed Parameter: Out5_Y0
   * Referenced by: '<S22>/Out5'
   */
  0.0,

  /* Expression: 0
   * Referenced by: '<S22>/Constant4'
   */
  0.0,

  /* Expression: zeros(16,1)
   * Referenced by: '<S99>/SwitchCurrents'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 },

  /* Computed Parameter: DiscreteTimeIntegrator_gainval
   * Referenced by: '<S90>/Discrete-Time Integrator'
   */
  5.0E-6,

  /* Expression: 0
   * Referenced by: '<S90>/Discrete-Time Integrator'
   */
  0.0,

  /* Expression: S.D
   * Referenced by: '<S97>/State-Space'
   */
  { -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0,
    50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0,
    -50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0,
    0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, 50000.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    50000.0, 0.0, -50000.0, 0.0, -50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0,
    50000.0, 0.0, 0.0, -50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0,
    0.0, 50000.0, 0.0, 50000.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0,
    -50000.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, 0.0,
    0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0,
    -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0,
    -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0,
    0.0, 0.0, 0.0, -50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0,
    0.0, -50000.0, 0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 50000.0,
    0.0, -50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 50000.0,
    0.0, 50000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 50000.0, 0.0, -50000.0, 0.0, 0.0, 50000.0, 0.0, 50000.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 0.0, 50000.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 50000.0, 0.0, -50000.0, 0.0, -50000.0, 0.0, -50000.0, 0.0, 0.0, 50000.0,
    -50000.0, -50000.0, 50000.0, 0.0, 0.0, 0.0, 0.0, -50000.0, 50000.0, 50000.0,
    -50000.0, 0.0, 0.0, 0.0, 0.0, -100000.0, 0.0, -100000.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 50000.0, -50000.0, -50000.0, 50000.0, 0.0, 0.0, 0.0, 0.0,
    -50000.0, 50000.0, 50000.0, -50000.0, 0.0, -100000.0, 0.0, -100000.0, 0.0,
    1.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, -0.5, -0.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5,
    0.5, 0.5, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0 },

  /* Expression: 1
   * Referenced by: '<S87>/eliminate warning with bus selector'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S1>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S2>/do not delete this gain'
   */
  1.0,

  /* Expression: 30
   * Referenced by: '<Root>/Constant8'
   */
  30.0,

  /* Computed Parameter: DiscreteTimeIntegrator1_gainval
   * Referenced by: '<S89>/Discrete-Time Integrator1'
   */
  5.0E-6,

  /* Expression: w0
   * Referenced by: '<S89>/Discrete-Time Integrator1'
   */
  0.0,

  /* Computed Parameter: DiscreteTimeIntegrator_gainva_k
   * Referenced by: '<S89>/Discrete-Time Integrator'
   */
  5.0E-6,

  /* Expression: SM.theta0
   * Referenced by: '<S89>/Discrete-Time Integrator'
   */
  0.0,

  /* Expression: 250/1000
   * Referenced by: '<Root>/Constant'
   */
  0.25,

  /* Expression: 4
   * Referenced by: '<Root>/Gain5'
   */
  4.0,

  /* Expression: 0
   * Referenced by: '<Root>/Constant9'
   */
  0.0,

  /* Expression: sps.Delay
   * Referenced by: '<S91>/Constant3'
   */
  0.0,

  /* Expression: sps.Period
   * Referenced by: '<S91>/Constant1'
   */
  5.0E-5,

  /* Expression: sps.Freq
   * Referenced by: '<S91>/1\ib1'
   */
  20000.0,

  /* Expression: [0 2 0]
   * Referenced by: '<S91>/1-D Lookup Table'
   */
  { 0.0, 2.0, 0.0 },

  /* Expression: [0 .5 1]
   * Referenced by: '<S91>/1-D Lookup Table'
   */
  { 0.0, 0.5, 1.0 },

  /* Expression: 1
   * Referenced by: '<S91>/Constant2'
   */
  1.0,

  /* Expression: 0
   * Referenced by: '<Root>/Constant11'
   */
  0.0,

  /* Expression: 24
   * Referenced by: '<Root>/Constant1'
   */
  24.0,

  /* Expression: SM.p
   * Referenced by: '<S88>/p'
   */
  50.0,

  /* Expression: [0  -pi/2]
   * Referenced by: '<S87>/Constant'
   */
  { 0.0, -1.5707963267948966 },

  /* Expression: 2*pi
   * Referenced by: '<S88>/Constant1'
   */
  6.2831853071795862,

  /* Expression: -SM.p*Psim
   * Referenced by: '<S88>/pxPsim'
   */
  -0.33,

  /* Expression: 4
   * Referenced by: '<S88>/2'
   */
  4.0,

  /* Expression: -Tdm
   * Referenced by: '<S90>/Tdm'
   */
  -0.0,

  /* Expression: B
   * Referenced by: '<S89>/B'
   */
  0.0005,

  /* Expression: 1/J
   * Referenced by: '<S89>/1//J'
   */
  212765.95744680852,

  /* Expression: R
   * Referenced by: '<S90>/R'
   */
  3.4,

  /* Expression: 1/L
   * Referenced by: '<S90>/1//L'
   */
  166.66666666666666,

  /* Expression: 1
   * Referenced by: '<S24>/do not delete this gain'
   */
  1.0,

  /* Expression: 1
   * Referenced by: '<S25>/do not delete this gain'
   */
  1.0
};
