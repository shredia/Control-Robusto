/*
 * Simulacion_stepper.h
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

#ifndef Simulacion_stepper_h_
#define Simulacion_stepper_h_
#ifndef Simulacion_stepper_COMMON_INCLUDES_
#define Simulacion_stepper_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "rt_logging.h"
#include "rt_nonfinite.h"
#include "math.h"
#endif                                 /* Simulacion_stepper_COMMON_INCLUDES_ */

#include "Simulacion_stepper_types.h"
#include "rtGetNaN.h"
#include "rtGetInf.h"
#include <float.h>
#include <string.h>
#include <stddef.h>

/* Macros for accessing real-time model data structure */
#ifndef rtmGetFinalTime
#define rtmGetFinalTime(rtm)           ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetRTWLogInfo
#define rtmGetRTWLogInfo(rtm)          ((rtm)->rtwLogInfo)
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   ((rtm)->Timing.taskTime0)
#endif

#ifndef rtmGetTFinal
#define rtmGetTFinal(rtm)              ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                (&(rtm)->Timing.taskTime0)
#endif

/* Block signals for system '<Root>/MATLAB Function' */
typedef struct {
  real_T y;                            /* '<Root>/MATLAB Function' */
} B_MATLABFunction_Simulacion_s_T;

/* Block signals for system '<Root>/MATLAB Function1' */
typedef struct {
  real_T y;                            /* '<Root>/MATLAB Function1' */
} B_MATLABFunction1_Simulacion__T;

/* Block signals (default storage) */
typedef struct {
  real_T IphA[2];                      /* '<S90>/Discrete-Time Integrator' */
  real_T StateSpace_o1[22];            /* '<S97>/State-Space' */
  real_T StateSpace_o2[16];            /* '<S97>/State-Space' */
  real_T VphV[2];              /* '<S87>/eliminate warning with bus selector' */
  real_T donotdeletethisgain;          /* '<S1>/do not delete this gain' */
  real_T donotdeletethisgain_d;        /* '<S2>/do not delete this gain' */
  real_T DiscreteTimeIntegrator1;      /* '<S89>/Discrete-Time Integrator1' */
  real_T DiscreteTimeIntegrator;       /* '<S89>/Discrete-Time Integrator' */
  real_T Gain2;                        /* '<Root>/Gain2' */
  real_T Gain5;                        /* '<Root>/Gain5' */
  real_T Sin;                          /* '<Root>/Sin' */
  real_T Gain4;                        /* '<Root>/Gain4' */
  real_T Sum3;                         /* '<Root>/Sum3' */
  real_T alpha;                        /* '<S18>/alpha' */
  real_T beta;                         /* '<S18>/beta' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T DigitalClock;                 /* '<S91>/Digital Clock' */
  real_T Add1;                         /* '<S91>/Add1' */
  real_T MathFunction;                 /* '<S91>/Math Function' */
  real_T uib1;                         /* '<S91>/1\ib1' */
  real_T uDLookupTable;                /* '<S91>/1-D Lookup Table' */
  real_T Add3;                         /* '<S91>/Add3' */
  real_T Gain1;                        /* '<Root>/Gain1' */
  real_T alpha_p;                      /* '<S19>/alpha' */
  real_T beta_n;                       /* '<S19>/beta' */
  real_T Gain3;                        /* '<Root>/Gain3' */
  real_T Product;                      /* '<Root>/Product' */
  real_T Product1;                     /* '<Root>/Product1' */
  real_T Sum13;                        /* '<Root>/Sum13' */
  real_T Sqrt;                         /* '<Root>/Sqrt' */
  real_T Product2;                     /* '<Root>/Product2' */
  real_T Product3;                     /* '<Root>/Product3' */
  real_T Sum14;                        /* '<Root>/Sum14' */
  real_T Sqrt1;                        /* '<Root>/Sqrt1' */
  real_T Sum2;                         /* '<Root>/Sum2' */
  real_T Sum4;                         /* '<Root>/Sum4' */
  real_T Sum8;                         /* '<Root>/Sum8' */
  real_T Gain7;                        /* '<Root>/Gain7' */
  real_T Gain6;                        /* '<Root>/Gain6' */
  real_T Gain8;                        /* '<Root>/Gain8' */
  real_T Sum1;                         /* '<Root>/Sum1' */
  real_T Sum5;                         /* '<Root>/Sum5' */
  real_T Sum9;                         /* '<Root>/Sum9' */
  real_T p;                            /* '<S88>/p' */
  real_T Sum1_k[2];                    /* '<S88>/Sum1' */
  real_T MathFunction_h[2];            /* '<S88>/Math Function' */
  real_T TrigonometricFunction[2];     /* '<S88>/Trigonometric Function' */
  real_T pxPsim[2];                    /* '<S88>/pxPsim' */
  real_T Product_c[2];                 /* '<S90>/Product' */
  real_T SumofElements;                /* '<S90>/Sum of Elements' */
  real_T u;                            /* '<S88>/2' */
  real_T MathFunction1;                /* '<S88>/Math Function1' */
  real_T TrigonometricFunction1;       /* '<S88>/Trigonometric Function1' */
  real_T Tdm;                          /* '<S90>/Tdm' */
  real_T TeNm;                         /* '<S90>/Sum2' */
  real_T DataTypeConversion;           /* '<S36>/Data Type Conversion' */
  real_T DataTypeConversion_b;         /* '<S43>/Data Type Conversion' */
  real_T DataTypeConversion_n;         /* '<S50>/Data Type Conversion' */
  real_T DataTypeConversion_p;         /* '<S57>/Data Type Conversion' */
  real_T DataTypeConversion_k;         /* '<S64>/Data Type Conversion' */
  real_T DataTypeConversion_ky;        /* '<S71>/Data Type Conversion' */
  real_T DataTypeConversion_d;         /* '<S78>/Data Type Conversion' */
  real_T DataTypeConversion_i;         /* '<S85>/Data Type Conversion' */
  real_T Product1_m[2];                /* '<S88>/Product1' */
  real_T B;                            /* '<S89>/B' */
  real_T Sum1_l;                       /* '<S89>/Sum1' */
  real_T uJ;                           /* '<S89>/1//J' */
  real_T R[2];                         /* '<S90>/R' */
  real_T Sum[2];                       /* '<S90>/Sum' */
  real_T IphA_e[2];                    /* '<S90>/1//L' */
  real_T donotdeletethisgain_j;        /* '<S24>/do not delete this gain' */
  real_T donotdeletethisgain_n;        /* '<S25>/do not delete this gain' */
  real_T CCCodeBlock_o1;               /* '<S22>/C//C++ Code Block' */
  real_T CCCodeBlock_o2;               /* '<S22>/C//C++ Code Block' */
  real_T CCCodeBlock_o3;               /* '<S22>/C//C++ Code Block' */
  real_T CCCodeBlock_o4;               /* '<S22>/C//C++ Code Block' */
  real_T CCCodeBlock_o5;               /* '<S22>/C//C++ Code Block' */
  real_T TmpSignalConversionAtSFunctionI[2];/* '<S23>/MATLAB Function2' */
  real_T TmpSignalConversionAtSFunctio_j[2];/* '<S23>/MATLAB Function2' */
  real_T x_corr[4];                    /* '<S23>/MATLAB Function2' */
  real_T TmpSignalConversionAtSFunctio_n[2];/* '<S23>/MATLAB Function1' */
  real_T TmpSignalConversionAtSFunctio_h[2];/* '<S23>/MATLAB Function1' */
  real_T x_corr_c[5];                  /* '<S23>/MATLAB Function1' */
  real_T TmpSignalConversionAtSFunctio_f[2];/* '<S23>/MATLAB Function' */
  real_T TmpSignalConversionAtSFunctio_o[2];/* '<S23>/MATLAB Function' */
  real_T x_corr_a[4];                  /* '<S23>/MATLAB Function' */
  int16_T RelationalOperator;          /* '<Root>/Relational Operator' */
  int16_T RelationalOperator4;         /* '<Root>/Relational Operator4' */
  int16_T Sum6;                        /* '<Root>/Sum6' */
  int16_T RelationalOperator1;         /* '<Root>/Relational Operator1' */
  int16_T RelationalOperator5;         /* '<Root>/Relational Operator5' */
  int16_T Sum7;                        /* '<Root>/Sum7' */
  boolean_T NOT;                       /* '<Root>/NOT' */
  boolean_T NOT4;                      /* '<Root>/NOT4' */
  boolean_T NOT1;                      /* '<Root>/NOT1' */
  boolean_T NOT5;                      /* '<Root>/NOT5' */
  B_MATLABFunction_Simulacion_s_T sf_MATLABFunction4;/* '<Root>/MATLAB Function4' */
  B_MATLABFunction_Simulacion_s_T sf_MATLABFunction3;/* '<Root>/MATLAB Function3' */
  B_MATLABFunction1_Simulacion__T sf_MATLABFunction2;/* '<Root>/MATLAB Function2' */
  B_MATLABFunction1_Simulacion__T sf_MATLABFunction1;/* '<Root>/MATLAB Function1' */
  B_MATLABFunction_Simulacion_s_T sf_MATLABFunction;/* '<Root>/MATLAB Function' */
} B_Simulacion_stepper_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteTimeIntegrator_DSTATE[2];/* '<S90>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator1_DSTATE;/* '<S89>/Discrete-Time Integrator1' */
  real_T DiscreteTimeIntegrator_DSTATE_o;/* '<S89>/Discrete-Time Integrator' */
  real_T P[16];                        /* '<S23>/MATLAB Function2' */
  real_T alpham[21];                   /* '<S23>/MATLAB Function2' */
  real_T alphac[21];                   /* '<S23>/MATLAB Function2' */
  real_T Qk[16];                       /* '<S23>/MATLAB Function2' */
  real_T Rk[4];                        /* '<S23>/MATLAB Function2' */
  real_T h;                            /* '<S23>/MATLAB Function2' */
  real_T x_hat[5];                     /* '<S23>/MATLAB Function1' */
  real_T P_p[25];                      /* '<S23>/MATLAB Function1' */
  real_T x_hat_l[4];                   /* '<S23>/MATLAB Function' */
  real_T P_n[16];                      /* '<S23>/MATLAB Function' */
  struct {
    void *AS;
    void *BS;
    void *CS;
    void *DS;
    void *DX_COL;
    void *BD_COL;
    void *TMP1;
    void *TMP2;
    void *XTMP;
    void *SWITCH_STATUS;
    void *SWITCH_STATUS_INIT;
    void *SW_CHG;
    void *G_STATE;
    void *USWLAST;
    void *XKM12;
    void *XKP12;
    void *XLAST;
    void *ULAST;
    void *IDX_SW_CHG;
    void *Y_SWITCH;
    void *SWITCH_TYPES;
    void *IDX_OUT_SW;
    void *SWITCH_TOPO_SAVED_IDX;
    void *SWITCH_MAP;
  } StateSpace_PWORK;                  /* '<S97>/State-Space' */

  uint32_T m_bpIndex;                  /* '<S91>/1-D Lookup Table' */
  int_T StateSpace_IWORK[11];          /* '<S97>/State-Space' */
} DW_Simulacion_stepper_T;

/* Parameters (default storage) */
struct P_Simulacion_stepper_T_ {
  real_T B;                            /* Variable: B
                                        * Referenced by:
                                        *   '<S22>/Constant'
                                        *   '<S23>/Constant'
                                        *   '<S23>/Constant25'
                                        */
  real_T J;                            /* Variable: J
                                        * Referenced by:
                                        *   '<S23>/Constant22'
                                        *   '<S23>/Constant29'
                                        */
  real_T Ki_corriente;                 /* Variable: Ki_corriente
                                        * Referenced by: '<S22>/Constant3'
                                        */
  real_T Ki_w;                         /* Variable: Ki_w
                                        * Referenced by: '<S22>/Constant11'
                                        */
  real_T Kp_corriente;                 /* Variable: Kp_corriente
                                        * Referenced by: '<S22>/Constant2'
                                        */
  real_T Kp_w;                         /* Variable: Kp_w
                                        * Referenced by: '<S22>/Constant1'
                                        */
  real_T Kt;                           /* Variable: Kt
                                        * Referenced by:
                                        *   '<S22>/Constant10'
                                        *   '<S23>/Constant20'
                                        *   '<S23>/Constant27'
                                        */
  real_T L;                            /* Variable: L
                                        * Referenced by:
                                        *   '<S22>/Constant8'
                                        *   '<S22>/Constant9'
                                        *   '<S23>/Constant23'
                                        *   '<S23>/Constant30'
                                        */
  real_T P;                            /* Variable: P
                                        * Referenced by:
                                        *   '<Root>/Gain2'
                                        *   '<Root>/Gain3'
                                        *   '<Root>/Gain6'
                                        *   '<Root>/Gain7'
                                        *   '<Root>/Gain8'
                                        *   '<S23>/Constant24'
                                        *   '<S23>/Constant31'
                                        */
  real_T R;                            /* Variable: R
                                        * Referenced by:
                                        *   '<S22>/Constant7'
                                        *   '<S23>/Constant21'
                                        *   '<S23>/Constant28'
                                        */
  real_T Tdm;                          /* Variable: Tdm
                                        * Referenced by: '<Root>/Gain4'
                                        */
  real_T Vdc;                          /* Variable: Vdc
                                        * Referenced by:
                                        *   '<Root>/Gain'
                                        *   '<Root>/Gain1'
                                        *   '<S22>/Constant6'
                                        *   '<S29>/DC'
                                        *   '<S30>/DC'
                                        */
  real_T sample_time;                  /* Variable: sample_time
                                        * Referenced by:
                                        *   '<S22>/Constant5'
                                        *   '<S23>/Constant19'
                                        *   '<S23>/Constant26'
                                        */
  real_T ekf_4var1_Y0;                 /* Computed Parameter: ekf_4var1_Y0
                                        * Referenced by: '<S23>/ekf_4var1'
                                        */
  real_T ekf_5var_Y0;                  /* Computed Parameter: ekf_5var_Y0
                                        * Referenced by: '<S23>/ekf_5var'
                                        */
  real_T Out1_Y0;                      /* Computed Parameter: Out1_Y0
                                        * Referenced by: '<S22>/Out1'
                                        */
  real_T Out2_Y0;                      /* Computed Parameter: Out2_Y0
                                        * Referenced by: '<S22>/Out2'
                                        */
  real_T Out3_Y0;                      /* Computed Parameter: Out3_Y0
                                        * Referenced by: '<S22>/Out3'
                                        */
  real_T Out4_Y0;                      /* Computed Parameter: Out4_Y0
                                        * Referenced by: '<S22>/Out4'
                                        */
  real_T Out5_Y0;                      /* Computed Parameter: Out5_Y0
                                        * Referenced by: '<S22>/Out5'
                                        */
  real_T Constant4_Value;              /* Expression: 0
                                        * Referenced by: '<S22>/Constant4'
                                        */
  real_T SwitchCurrents_Value[16];     /* Expression: zeros(16,1)
                                        * Referenced by: '<S99>/SwitchCurrents'
                                        */
  real_T DiscreteTimeIntegrator_gainval;
                           /* Computed Parameter: DiscreteTimeIntegrator_gainval
                            * Referenced by: '<S90>/Discrete-Time Integrator'
                            */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: 0
                                        * Referenced by: '<S90>/Discrete-Time Integrator'
                                        */
  real_T StateSpace_DS_param[440];     /* Expression: S.D
                                        * Referenced by: '<S97>/State-Space'
                                        */
  real_T eliminatewarningwithbusselector;/* Expression: 1
                                          * Referenced by: '<S87>/eliminate warning with bus selector'
                                          */
  real_T donotdeletethisgain_Gain;     /* Expression: 1
                                        * Referenced by: '<S1>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_l;   /* Expression: 1
                                        * Referenced by: '<S2>/do not delete this gain'
                                        */
  real_T Constant8_Value;              /* Expression: 30
                                        * Referenced by: '<Root>/Constant8'
                                        */
  real_T DiscreteTimeIntegrator1_gainval;
                          /* Computed Parameter: DiscreteTimeIntegrator1_gainval
                           * Referenced by: '<S89>/Discrete-Time Integrator1'
                           */
  real_T DiscreteTimeIntegrator1_IC;   /* Expression: w0
                                        * Referenced by: '<S89>/Discrete-Time Integrator1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_k;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_k
                           * Referenced by: '<S89>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_l;  /* Expression: SM.theta0
                                        * Referenced by: '<S89>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value;               /* Expression: 250/1000
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Gain5_Gain;                   /* Expression: 4
                                        * Referenced by: '<Root>/Gain5'
                                        */
  real_T Constant9_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant9'
                                        */
  real_T Constant3_Value;              /* Expression: sps.Delay
                                        * Referenced by: '<S91>/Constant3'
                                        */
  real_T Constant1_Value;              /* Expression: sps.Period
                                        * Referenced by: '<S91>/Constant1'
                                        */
  real_T uib1_Gain;                    /* Expression: sps.Freq
                                        * Referenced by: '<S91>/1\ib1'
                                        */
  real_T uDLookupTable_tableData[3];   /* Expression: [0 2 0]
                                        * Referenced by: '<S91>/1-D Lookup Table'
                                        */
  real_T uDLookupTable_bp01Data[3];    /* Expression: [0 .5 1]
                                        * Referenced by: '<S91>/1-D Lookup Table'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S91>/Constant2'
                                        */
  real_T Constant11_Value;             /* Expression: 0
                                        * Referenced by: '<Root>/Constant11'
                                        */
  real_T Constant1_Value_m;            /* Expression: 24
                                        * Referenced by: '<Root>/Constant1'
                                        */
  real_T p_Gain;                       /* Expression: SM.p
                                        * Referenced by: '<S88>/p'
                                        */
  real_T Constant_Value_f[2];          /* Expression: [0  -pi/2]
                                        * Referenced by: '<S87>/Constant'
                                        */
  real_T Constant1_Value_md;           /* Expression: 2*pi
                                        * Referenced by: '<S88>/Constant1'
                                        */
  real_T pxPsim_Gain;                  /* Expression: -SM.p*Psim
                                        * Referenced by: '<S88>/pxPsim'
                                        */
  real_T u_Gain;                       /* Expression: 4
                                        * Referenced by: '<S88>/2'
                                        */
  real_T Tdm_Gain;                     /* Expression: -Tdm
                                        * Referenced by: '<S90>/Tdm'
                                        */
  real_T B_Gain;                       /* Expression: B
                                        * Referenced by: '<S89>/B'
                                        */
  real_T uJ_Gain;                      /* Expression: 1/J
                                        * Referenced by: '<S89>/1//J'
                                        */
  real_T R_Gain;                       /* Expression: R
                                        * Referenced by: '<S90>/R'
                                        */
  real_T uL_Gain;                      /* Expression: 1/L
                                        * Referenced by: '<S90>/1//L'
                                        */
  real_T donotdeletethisgain_Gain_j;   /* Expression: 1
                                        * Referenced by: '<S24>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_i;   /* Expression: 1
                                        * Referenced by: '<S25>/do not delete this gain'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Simulacion_stepper_T {
  const char_T *errorStatus;
  RTWLogInfo *rtwLogInfo;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T taskTime0;
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    struct {
      uint8_T TID[2];
    } TaskCounters;

    time_T tFinal;
    boolean_T stopRequestedFlag;
  } Timing;
};

/* Block parameters (default storage) */
extern P_Simulacion_stepper_T Simulacion_stepper_P;

/* Block signals (default storage) */
extern B_Simulacion_stepper_T Simulacion_stepper_B;

/* Block states (default storage) */
extern DW_Simulacion_stepper_T Simulacion_stepper_DW;

/* Model entry point functions */
extern void Simulacion_stepper_initialize(void);
extern void Simulacion_stepper_step(void);
extern void Simulacion_stepper_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Simulacion_stepper_T *const Simulacion_stepper_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S10>/Add' : Unused code path elimination
 * Block '<S34>/0 1' : Unused code path elimination
 * Block '<S34>/Gain' : Unused code path elimination
 * Block '<S34>/Saturation' : Unused code path elimination
 * Block '<S34>/Sum' : Unused code path elimination
 * Block '<S34>/Switch' : Unused code path elimination
 * Block '<S34>/eee' : Unused code path elimination
 * Block '<S36>/0 1' : Unused code path elimination
 * Block '<S36>/1//Rsw' : Unused code path elimination
 * Block '<S36>/Switch' : Unused code path elimination
 * Block '<S11>/Add' : Unused code path elimination
 * Block '<S41>/0 1' : Unused code path elimination
 * Block '<S41>/Gain' : Unused code path elimination
 * Block '<S41>/Saturation' : Unused code path elimination
 * Block '<S41>/Sum' : Unused code path elimination
 * Block '<S41>/Switch' : Unused code path elimination
 * Block '<S41>/eee' : Unused code path elimination
 * Block '<S43>/0 1' : Unused code path elimination
 * Block '<S43>/1//Rsw' : Unused code path elimination
 * Block '<S43>/Switch' : Unused code path elimination
 * Block '<S12>/Add' : Unused code path elimination
 * Block '<S48>/0 1' : Unused code path elimination
 * Block '<S48>/Gain' : Unused code path elimination
 * Block '<S48>/Saturation' : Unused code path elimination
 * Block '<S48>/Sum' : Unused code path elimination
 * Block '<S48>/Switch' : Unused code path elimination
 * Block '<S48>/eee' : Unused code path elimination
 * Block '<S50>/0 1' : Unused code path elimination
 * Block '<S50>/1//Rsw' : Unused code path elimination
 * Block '<S50>/Switch' : Unused code path elimination
 * Block '<S13>/Add' : Unused code path elimination
 * Block '<S55>/0 1' : Unused code path elimination
 * Block '<S55>/Gain' : Unused code path elimination
 * Block '<S55>/Saturation' : Unused code path elimination
 * Block '<S55>/Sum' : Unused code path elimination
 * Block '<S55>/Switch' : Unused code path elimination
 * Block '<S55>/eee' : Unused code path elimination
 * Block '<S57>/0 1' : Unused code path elimination
 * Block '<S57>/1//Rsw' : Unused code path elimination
 * Block '<S57>/Switch' : Unused code path elimination
 * Block '<S14>/Add' : Unused code path elimination
 * Block '<S62>/0 1' : Unused code path elimination
 * Block '<S62>/Gain' : Unused code path elimination
 * Block '<S62>/Saturation' : Unused code path elimination
 * Block '<S62>/Sum' : Unused code path elimination
 * Block '<S62>/Switch' : Unused code path elimination
 * Block '<S62>/eee' : Unused code path elimination
 * Block '<S64>/0 1' : Unused code path elimination
 * Block '<S64>/1//Rsw' : Unused code path elimination
 * Block '<S64>/Switch' : Unused code path elimination
 * Block '<S15>/Add' : Unused code path elimination
 * Block '<S69>/0 1' : Unused code path elimination
 * Block '<S69>/Gain' : Unused code path elimination
 * Block '<S69>/Saturation' : Unused code path elimination
 * Block '<S69>/Sum' : Unused code path elimination
 * Block '<S69>/Switch' : Unused code path elimination
 * Block '<S69>/eee' : Unused code path elimination
 * Block '<S71>/0 1' : Unused code path elimination
 * Block '<S71>/1//Rsw' : Unused code path elimination
 * Block '<S71>/Switch' : Unused code path elimination
 * Block '<S16>/Add' : Unused code path elimination
 * Block '<S76>/0 1' : Unused code path elimination
 * Block '<S76>/Gain' : Unused code path elimination
 * Block '<S76>/Saturation' : Unused code path elimination
 * Block '<S76>/Sum' : Unused code path elimination
 * Block '<S76>/Switch' : Unused code path elimination
 * Block '<S76>/eee' : Unused code path elimination
 * Block '<S78>/0 1' : Unused code path elimination
 * Block '<S78>/1//Rsw' : Unused code path elimination
 * Block '<S78>/Switch' : Unused code path elimination
 * Block '<S17>/Add' : Unused code path elimination
 * Block '<S83>/0 1' : Unused code path elimination
 * Block '<S83>/Gain' : Unused code path elimination
 * Block '<S83>/Saturation' : Unused code path elimination
 * Block '<S83>/Sum' : Unused code path elimination
 * Block '<S83>/Switch' : Unused code path elimination
 * Block '<S83>/eee' : Unused code path elimination
 * Block '<S85>/0 1' : Unused code path elimination
 * Block '<S85>/1//Rsw' : Unused code path elimination
 * Block '<S85>/Switch' : Unused code path elimination
 * Block '<S18>/0' : Unused code path elimination
 * Block '<S19>/0' : Unused code path elimination
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Simulacion_stepper'
 * '<S1>'   : 'Simulacion_stepper/Current Measurement'
 * '<S2>'   : 'Simulacion_stepper/Current Measurement1'
 * '<S3>'   : 'Simulacion_stepper/DC Voltage Source'
 * '<S4>'   : 'Simulacion_stepper/DC Voltage Source1'
 * '<S5>'   : 'Simulacion_stepper/MATLAB Function'
 * '<S6>'   : 'Simulacion_stepper/MATLAB Function1'
 * '<S7>'   : 'Simulacion_stepper/MATLAB Function2'
 * '<S8>'   : 'Simulacion_stepper/MATLAB Function3'
 * '<S9>'   : 'Simulacion_stepper/MATLAB Function4'
 * '<S10>'  : 'Simulacion_stepper/Mosfet'
 * '<S11>'  : 'Simulacion_stepper/Mosfet1'
 * '<S12>'  : 'Simulacion_stepper/Mosfet2'
 * '<S13>'  : 'Simulacion_stepper/Mosfet3'
 * '<S14>'  : 'Simulacion_stepper/Mosfet4'
 * '<S15>'  : 'Simulacion_stepper/Mosfet5'
 * '<S16>'  : 'Simulacion_stepper/Mosfet6'
 * '<S17>'  : 'Simulacion_stepper/Mosfet7'
 * '<S18>'  : 'Simulacion_stepper/Park to Clarke Angle Transform'
 * '<S19>'  : 'Simulacion_stepper/Park to Clarke Angle Transform1'
 * '<S20>'  : 'Simulacion_stepper/Stepper Motor'
 * '<S21>'  : 'Simulacion_stepper/Triangle Generator'
 * '<S22>'  : 'Simulacion_stepper/Triggered Subsystem'
 * '<S23>'  : 'Simulacion_stepper/Triggered Subsystem1'
 * '<S24>'  : 'Simulacion_stepper/Voltage Measurement'
 * '<S25>'  : 'Simulacion_stepper/Voltage Measurement1'
 * '<S26>'  : 'Simulacion_stepper/powergui'
 * '<S27>'  : 'Simulacion_stepper/Current Measurement/Model'
 * '<S28>'  : 'Simulacion_stepper/Current Measurement1/Model'
 * '<S29>'  : 'Simulacion_stepper/DC Voltage Source/Model'
 * '<S30>'  : 'Simulacion_stepper/DC Voltage Source1/Model'
 * '<S31>'  : 'Simulacion_stepper/Mosfet/Diode'
 * '<S32>'  : 'Simulacion_stepper/Mosfet/Ideal Switch'
 * '<S33>'  : 'Simulacion_stepper/Mosfet/Measurement list'
 * '<S34>'  : 'Simulacion_stepper/Mosfet/Diode/Model'
 * '<S35>'  : 'Simulacion_stepper/Mosfet/Diode/Model/Measurement list'
 * '<S36>'  : 'Simulacion_stepper/Mosfet/Ideal Switch/Model'
 * '<S37>'  : 'Simulacion_stepper/Mosfet/Ideal Switch/Model/Measurement list'
 * '<S38>'  : 'Simulacion_stepper/Mosfet1/Diode'
 * '<S39>'  : 'Simulacion_stepper/Mosfet1/Ideal Switch'
 * '<S40>'  : 'Simulacion_stepper/Mosfet1/Measurement list'
 * '<S41>'  : 'Simulacion_stepper/Mosfet1/Diode/Model'
 * '<S42>'  : 'Simulacion_stepper/Mosfet1/Diode/Model/Measurement list'
 * '<S43>'  : 'Simulacion_stepper/Mosfet1/Ideal Switch/Model'
 * '<S44>'  : 'Simulacion_stepper/Mosfet1/Ideal Switch/Model/Measurement list'
 * '<S45>'  : 'Simulacion_stepper/Mosfet2/Diode'
 * '<S46>'  : 'Simulacion_stepper/Mosfet2/Ideal Switch'
 * '<S47>'  : 'Simulacion_stepper/Mosfet2/Measurement list'
 * '<S48>'  : 'Simulacion_stepper/Mosfet2/Diode/Model'
 * '<S49>'  : 'Simulacion_stepper/Mosfet2/Diode/Model/Measurement list'
 * '<S50>'  : 'Simulacion_stepper/Mosfet2/Ideal Switch/Model'
 * '<S51>'  : 'Simulacion_stepper/Mosfet2/Ideal Switch/Model/Measurement list'
 * '<S52>'  : 'Simulacion_stepper/Mosfet3/Diode'
 * '<S53>'  : 'Simulacion_stepper/Mosfet3/Ideal Switch'
 * '<S54>'  : 'Simulacion_stepper/Mosfet3/Measurement list'
 * '<S55>'  : 'Simulacion_stepper/Mosfet3/Diode/Model'
 * '<S56>'  : 'Simulacion_stepper/Mosfet3/Diode/Model/Measurement list'
 * '<S57>'  : 'Simulacion_stepper/Mosfet3/Ideal Switch/Model'
 * '<S58>'  : 'Simulacion_stepper/Mosfet3/Ideal Switch/Model/Measurement list'
 * '<S59>'  : 'Simulacion_stepper/Mosfet4/Diode'
 * '<S60>'  : 'Simulacion_stepper/Mosfet4/Ideal Switch'
 * '<S61>'  : 'Simulacion_stepper/Mosfet4/Measurement list'
 * '<S62>'  : 'Simulacion_stepper/Mosfet4/Diode/Model'
 * '<S63>'  : 'Simulacion_stepper/Mosfet4/Diode/Model/Measurement list'
 * '<S64>'  : 'Simulacion_stepper/Mosfet4/Ideal Switch/Model'
 * '<S65>'  : 'Simulacion_stepper/Mosfet4/Ideal Switch/Model/Measurement list'
 * '<S66>'  : 'Simulacion_stepper/Mosfet5/Diode'
 * '<S67>'  : 'Simulacion_stepper/Mosfet5/Ideal Switch'
 * '<S68>'  : 'Simulacion_stepper/Mosfet5/Measurement list'
 * '<S69>'  : 'Simulacion_stepper/Mosfet5/Diode/Model'
 * '<S70>'  : 'Simulacion_stepper/Mosfet5/Diode/Model/Measurement list'
 * '<S71>'  : 'Simulacion_stepper/Mosfet5/Ideal Switch/Model'
 * '<S72>'  : 'Simulacion_stepper/Mosfet5/Ideal Switch/Model/Measurement list'
 * '<S73>'  : 'Simulacion_stepper/Mosfet6/Diode'
 * '<S74>'  : 'Simulacion_stepper/Mosfet6/Ideal Switch'
 * '<S75>'  : 'Simulacion_stepper/Mosfet6/Measurement list'
 * '<S76>'  : 'Simulacion_stepper/Mosfet6/Diode/Model'
 * '<S77>'  : 'Simulacion_stepper/Mosfet6/Diode/Model/Measurement list'
 * '<S78>'  : 'Simulacion_stepper/Mosfet6/Ideal Switch/Model'
 * '<S79>'  : 'Simulacion_stepper/Mosfet6/Ideal Switch/Model/Measurement list'
 * '<S80>'  : 'Simulacion_stepper/Mosfet7/Diode'
 * '<S81>'  : 'Simulacion_stepper/Mosfet7/Ideal Switch'
 * '<S82>'  : 'Simulacion_stepper/Mosfet7/Measurement list'
 * '<S83>'  : 'Simulacion_stepper/Mosfet7/Diode/Model'
 * '<S84>'  : 'Simulacion_stepper/Mosfet7/Diode/Model/Measurement list'
 * '<S85>'  : 'Simulacion_stepper/Mosfet7/Ideal Switch/Model'
 * '<S86>'  : 'Simulacion_stepper/Mosfet7/Ideal Switch/Model/Measurement list'
 * '<S87>'  : 'Simulacion_stepper/Stepper Motor/Model'
 * '<S88>'  : 'Simulacion_stepper/Stepper Motor/Model/FEM'
 * '<S89>'  : 'Simulacion_stepper/Stepper Motor/Model/Mechanical'
 * '<S90>'  : 'Simulacion_stepper/Stepper Motor/Model/windings'
 * '<S91>'  : 'Simulacion_stepper/Triangle Generator/Model'
 * '<S92>'  : 'Simulacion_stepper/Triggered Subsystem1/MATLAB Function'
 * '<S93>'  : 'Simulacion_stepper/Triggered Subsystem1/MATLAB Function1'
 * '<S94>'  : 'Simulacion_stepper/Triggered Subsystem1/MATLAB Function2'
 * '<S95>'  : 'Simulacion_stepper/Voltage Measurement/Model'
 * '<S96>'  : 'Simulacion_stepper/Voltage Measurement1/Model'
 * '<S97>'  : 'Simulacion_stepper/powergui/EquivalentModel1'
 * '<S98>'  : 'Simulacion_stepper/powergui/EquivalentModel1/Gates'
 * '<S99>'  : 'Simulacion_stepper/powergui/EquivalentModel1/Sources'
 * '<S100>' : 'Simulacion_stepper/powergui/EquivalentModel1/Status'
 * '<S101>' : 'Simulacion_stepper/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* Simulacion_stepper_h_ */
