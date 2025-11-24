/*
 * Simulacion_stepper_private.h
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

#ifndef Simulacion_stepper_private_h_
#define Simulacion_stepper_private_h_
#include "rtwtypes.h"
#include "builtin_typeid_types.h"
#include "multiword_types.h"
#include "Simulacion_stepper.h"
#include "Simulacion_stepper_types.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include <math.h>
#include <stdlib.h>

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetTFinal
#define rtmSetTFinal(rtm, val)         ((rtm)->Timing.tFinal = (val))
#endif

#ifndef CodeFormat
#define CodeFormat                     S-Function
#else
#undef CodeFormat
#define CodeFormat                     S-Function
#endif

#ifndef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#else
#undef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#endif

#ifndef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#else
#undef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#endif

#ifndef RTW_GENERATED_S_FUNCTION
#define RTW_GENERATED_S_FUNCTION
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        NULL
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)
#endif

#if !defined(RTW_SFUNCTION_DEFINES)
#define RTW_SFUNCTION_DEFINES
#ifndef _RTW_COMMON_DEFINES_
#define _RTW_COMMON_DEFINES_
#endif
#endif

#ifdef __cplusplus
#define SFB_EXTERN_C                   extern "C"
#else
#define SFB_EXTERN_C                   extern
#endif

SFB_EXTERN_C void PI_dq_Start_wrapper(void);
SFB_EXTERN_C void PI_dq_Outputs_wrapper(const real_T *Ia,
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
  const real_T *Tl_Tdm,
  const real_T *B,
  real_T *Ud,
  real_T *Uq,
  real_T *Id,
  real_T *Iq,
  real_T *Iq_ref_aux);
SFB_EXTERN_C void PI_dq_Terminate_wrapper(void);

#undef SFB_EXTERN_C

extern real_T rt_remd_snf(real_T u0, real_T u1);
extern real_T rt_modd_snf(real_T u0, real_T u1);
extern real_T look1_pbinlxmpw(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T prevIndex[], uint32_T maxIndex);
extern int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator);
extern void Simulacion_stepp_MATLABFunction(real_T rtu_u,
  B_MATLABFunction_Simulacion_s_T *localB);
extern void Simulacion_step_MATLABFunction1(real_T rtu_u,
  B_MATLABFunction1_Simulacion__T *localB);

#endif                                 /* Simulacion_stepper_private_h_ */
