/*
 * Simulacion_stepper.h
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

/* Block signals (default storage) */
typedef struct {
  real_T IphA[2];                      /* '<S85>/Discrete-Time Integrator' */
  real_T StateSpace_o1[22];            /* '<S89>/State-Space' */
  real_T StateSpace_o2[16];            /* '<S89>/State-Space' */
  real_T VphV[2];              /* '<S82>/eliminate warning with bus selector' */
  real_T donotdeletethisgain;          /* '<S1>/do not delete this gain' */
  real_T donotdeletethisgain_d;        /* '<S2>/do not delete this gain' */
  real_T DiscreteTimeIntegrator1;      /* '<S84>/Discrete-Time Integrator1' */
  real_T DiscreteTimeIntegrator;       /* '<S84>/Discrete-Time Integrator' */
  real_T Gain2;                        /* '<Root>/Gain2' */
  real_T Gain5;                        /* '<Root>/Gain5' */
  real_T Sin;                          /* '<Root>/Sin' */
  real_T Gain4;                        /* '<Root>/Gain4' */
  real_T Sum3;                         /* '<Root>/Sum3' */
  real_T alpha;                        /* '<S13>/alpha' */
  real_T beta;                         /* '<S13>/beta' */
  real_T Gain;                         /* '<Root>/Gain' */
  real_T DigitalClock;                 /* '<S86>/Digital Clock' */
  real_T Add1;                         /* '<S86>/Add1' */
  real_T MathFunction;                 /* '<S86>/Math Function' */
  real_T uib1;                         /* '<S86>/1\ib1' */
  real_T uDLookupTable;                /* '<S86>/1-D Lookup Table' */
  real_T Add3;                         /* '<S86>/Add3' */
  real_T Gain1;                        /* '<Root>/Gain1' */
  real_T alpha_p;                      /* '<S14>/alpha' */
  real_T beta_n;                       /* '<S14>/beta' */
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
  real_T Sum1;                         /* '<Root>/Sum1' */
  real_T Sum5;                         /* '<Root>/Sum5' */
  real_T Sum9;                         /* '<Root>/Sum9' */
  real_T Mod;                          /* '<Root>/Mod' */
  real_T Mod2;                         /* '<Root>/Mod2' */
  real_T Mod3;                         /* '<Root>/Mod3' */
  real_T Mod1;                         /* '<Root>/Mod1' */
  real_T Gain9;                        /* '<Root>/Gain9' */
  real_T p;                            /* '<S83>/p' */
  real_T Sum1_k[2];                    /* '<S83>/Sum1' */
  real_T MathFunction_h[2];            /* '<S83>/Math Function' */
  real_T TrigonometricFunction[2];     /* '<S83>/Trigonometric Function' */
  real_T pxPsim[2];                    /* '<S83>/pxPsim' */
  real_T Product_c[2];                 /* '<S85>/Product' */
  real_T SumofElements;                /* '<S85>/Sum of Elements' */
  real_T u;                            /* '<S83>/2' */
  real_T MathFunction1;                /* '<S83>/Math Function1' */
  real_T TrigonometricFunction1;       /* '<S83>/Trigonometric Function1' */
  real_T Tdm;                          /* '<S85>/Tdm' */
  real_T TeNm;                         /* '<S85>/Sum2' */
  real_T DataTypeConversion;           /* '<S31>/Data Type Conversion' */
  real_T DataTypeConversion_b;         /* '<S38>/Data Type Conversion' */
  real_T DataTypeConversion_n;         /* '<S45>/Data Type Conversion' */
  real_T DataTypeConversion_p;         /* '<S52>/Data Type Conversion' */
  real_T DataTypeConversion_k;         /* '<S59>/Data Type Conversion' */
  real_T DataTypeConversion_ky;        /* '<S66>/Data Type Conversion' */
  real_T DataTypeConversion_d;         /* '<S73>/Data Type Conversion' */
  real_T DataTypeConversion_i;         /* '<S80>/Data Type Conversion' */
  real_T Product1_m[2];                /* '<S83>/Product1' */
  real_T B;                            /* '<S84>/B' */
  real_T Sum1_l;                       /* '<S84>/Sum1' */
  real_T uJ;                           /* '<S84>/1//J' */
  real_T R[2];                         /* '<S85>/R' */
  real_T Sum[2];                       /* '<S85>/Sum' */
  real_T IphA_e[2];                    /* '<S85>/1//L' */
  real_T donotdeletethisgain_j;        /* '<S19>/do not delete this gain' */
  real_T donotdeletethisgain_n;        /* '<S20>/do not delete this gain' */
  real_T CCCodeBlock_o1;               /* '<S17>/C//C++ Code Block' */
  real_T CCCodeBlock_o2;               /* '<S17>/C//C++ Code Block' */
  real_T CCCodeBlock_o3;               /* '<S17>/C//C++ Code Block' */
  real_T CCCodeBlock_o4;               /* '<S17>/C//C++ Code Block' */
  real_T CCCodeBlock_o5;               /* '<S17>/C//C++ Code Block' */
  real_T CCCodeBlock1[5];              /* '<S18>/C//C++ Code Block1' */
  real_T CCCodeBlock2[4];              /* '<S18>/C//C++ Code Block2' */
  real_T CCCodeBlock3_o1[4];           /* '<S18>/C//C++ Code Block3' */
  real_T CCCodeBlock3_o2[4];           /* '<S18>/C//C++ Code Block3' */
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
} B_Simulacion_stepper_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T DiscreteTimeIntegrator_DSTATE[2];/* '<S85>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator1_DSTATE;/* '<S84>/Discrete-Time Integrator1' */
  real_T DiscreteTimeIntegrator_DSTATE_o;/* '<S84>/Discrete-Time Integrator' */
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
  } StateSpace_PWORK;                  /* '<S89>/State-Space' */

  uint32_T m_bpIndex;                  /* '<S86>/1-D Lookup Table' */
  int_T StateSpace_IWORK[11];          /* '<S89>/State-Space' */
} DW_Simulacion_stepper_T;

/* Parameters (default storage) */
struct P_Simulacion_stepper_T_ {
  real_T B;                            /* Variable: B
                                        * Referenced by: '<S17>/Constant'
                                        */
  real_T J;                            /* Variable: J
                                        * Referenced by:
                                        *   '<S18>/Constant15'
                                        *   '<S18>/Constant3'
                                        *   '<S18>/Constant9'
                                        */
  real_T Ki_corriente;                 /* Variable: Ki_corriente
                                        * Referenced by: '<S17>/Constant3'
                                        */
  real_T Ki_w;                         /* Variable: Ki_w
                                        * Referenced by: '<S17>/Constant11'
                                        */
  real_T Kp_corriente;                 /* Variable: Kp_corriente
                                        * Referenced by: '<S17>/Constant2'
                                        */
  real_T Kp_w;                         /* Variable: Kp_w
                                        * Referenced by: '<S17>/Constant1'
                                        */
  real_T Kt;                           /* Variable: Kt
                                        * Referenced by:
                                        *   '<S17>/Constant10'
                                        *   '<S18>/Constant14'
                                        *   '<S18>/Constant2'
                                        *   '<S18>/Constant8'
                                        */
  real_T L;                            /* Variable: L
                                        * Referenced by:
                                        *   '<S17>/Constant8'
                                        *   '<S17>/Constant9'
                                        *   '<S18>/Constant1'
                                        *   '<S18>/Constant13'
                                        *   '<S18>/Constant7'
                                        */
  real_T P;                            /* Variable: P
                                        * Referenced by:
                                        *   '<Root>/Gain2'
                                        *   '<Root>/Gain3'
                                        *   '<Root>/Gain6'
                                        *   '<Root>/Gain7'
                                        *   '<Root>/Gain9'
                                        *   '<S18>/Constant11'
                                        *   '<S18>/Constant17'
                                        *   '<S18>/Constant5'
                                        */
  real_T R;                            /* Variable: R
                                        * Referenced by:
                                        *   '<S17>/Constant7'
                                        *   '<S18>/Constant12'
                                        *   '<S18>/Constant18'
                                        *   '<S18>/Constant6'
                                        */
  real_T Tdm;                          /* Variable: Tdm
                                        * Referenced by: '<Root>/Gain4'
                                        */
  real_T Vdc;                          /* Variable: Vdc
                                        * Referenced by:
                                        *   '<Root>/Gain'
                                        *   '<Root>/Gain1'
                                        *   '<S17>/Constant6'
                                        *   '<S24>/DC'
                                        *   '<S25>/DC'
                                        */
  real_T sample_time;                  /* Variable: sample_time
                                        * Referenced by:
                                        *   '<S17>/Constant5'
                                        *   '<S18>/Constant10'
                                        *   '<S18>/Constant16'
                                        *   '<S18>/Constant4'
                                        */
  real_T X_SPKF_Y0;                    /* Computed Parameter: X_SPKF_Y0
                                        * Referenced by: '<S18>/X_SPKF'
                                        */
  real_T ekf_5var_Y0;                  /* Computed Parameter: ekf_5var_Y0
                                        * Referenced by: '<S18>/ekf_5var'
                                        */
  real_T ekf_4var_Y0;                  /* Computed Parameter: ekf_4var_Y0
                                        * Referenced by: '<S18>/ekf_4var'
                                        */
  real_T Out1_Y0;                      /* Computed Parameter: Out1_Y0
                                        * Referenced by: '<S17>/Out1'
                                        */
  real_T Out2_Y0;                      /* Computed Parameter: Out2_Y0
                                        * Referenced by: '<S17>/Out2'
                                        */
  real_T Out3_Y0;                      /* Computed Parameter: Out3_Y0
                                        * Referenced by: '<S17>/Out3'
                                        */
  real_T Out4_Y0;                      /* Computed Parameter: Out4_Y0
                                        * Referenced by: '<S17>/Out4'
                                        */
  real_T Out5_Y0;                      /* Computed Parameter: Out5_Y0
                                        * Referenced by: '<S17>/Out5'
                                        */
  real_T Constant4_Value;              /* Expression: 0
                                        * Referenced by: '<S17>/Constant4'
                                        */
  real_T SwitchCurrents_Value[16];     /* Expression: zeros(16,1)
                                        * Referenced by: '<S91>/SwitchCurrents'
                                        */
  real_T DiscreteTimeIntegrator_gainval;
                           /* Computed Parameter: DiscreteTimeIntegrator_gainval
                            * Referenced by: '<S85>/Discrete-Time Integrator'
                            */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: 0
                                        * Referenced by: '<S85>/Discrete-Time Integrator'
                                        */
  real_T StateSpace_DS_param[440];     /* Expression: S.D
                                        * Referenced by: '<S89>/State-Space'
                                        */
  real_T eliminatewarningwithbusselector;/* Expression: 1
                                          * Referenced by: '<S82>/eliminate warning with bus selector'
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
                           * Referenced by: '<S84>/Discrete-Time Integrator1'
                           */
  real_T DiscreteTimeIntegrator1_IC;   /* Expression: w0
                                        * Referenced by: '<S84>/Discrete-Time Integrator1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_k;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_k
                           * Referenced by: '<S84>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_l;  /* Expression: SM.theta0
                                        * Referenced by: '<S84>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Gain5_Gain;                   /* Expression: 4
                                        * Referenced by: '<Root>/Gain5'
                                        */
  real_T Constant9_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/Constant9'
                                        */
  real_T Constant3_Value;              /* Expression: sps.Delay
                                        * Referenced by: '<S86>/Constant3'
                                        */
  real_T Constant1_Value;              /* Expression: sps.Period
                                        * Referenced by: '<S86>/Constant1'
                                        */
  real_T uib1_Gain;                    /* Expression: sps.Freq
                                        * Referenced by: '<S86>/1\ib1'
                                        */
  real_T uDLookupTable_tableData[3];   /* Expression: [0 2 0]
                                        * Referenced by: '<S86>/1-D Lookup Table'
                                        */
  real_T uDLookupTable_bp01Data[3];    /* Expression: [0 .5 1]
                                        * Referenced by: '<S86>/1-D Lookup Table'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S86>/Constant2'
                                        */
  real_T Constant11_Value;             /* Expression: 0
                                        * Referenced by: '<Root>/Constant11'
                                        */
  real_T Constant6_Value;              /* Expression: 2*pi
                                        * Referenced by: '<Root>/Constant6'
                                        */
  real_T Constant1_Value_j;            /* Expression: 2*pi
                                        * Referenced by: '<Root>/Constant1'
                                        */
  real_T Constant2_Value_m;            /* Expression: 2*pi
                                        * Referenced by: '<Root>/Constant2'
                                        */
  real_T Constant3_Value_j;            /* Expression: 2*pi
                                        * Referenced by: '<Root>/Constant3'
                                        */
  real_T p_Gain;                       /* Expression: SM.p
                                        * Referenced by: '<S83>/p'
                                        */
  real_T Constant_Value_f[2];          /* Expression: [0  -pi/2]
                                        * Referenced by: '<S82>/Constant'
                                        */
  real_T Constant1_Value_m;            /* Expression: 2*pi
                                        * Referenced by: '<S83>/Constant1'
                                        */
  real_T pxPsim_Gain;                  /* Expression: -SM.p*Psim
                                        * Referenced by: '<S83>/pxPsim'
                                        */
  real_T u_Gain;                       /* Expression: 4
                                        * Referenced by: '<S83>/2'
                                        */
  real_T Tdm_Gain;                     /* Expression: -Tdm
                                        * Referenced by: '<S85>/Tdm'
                                        */
  real_T B_Gain;                       /* Expression: B
                                        * Referenced by: '<S84>/B'
                                        */
  real_T uJ_Gain;                      /* Expression: 1/J
                                        * Referenced by: '<S84>/1//J'
                                        */
  real_T R_Gain;                       /* Expression: R
                                        * Referenced by: '<S85>/R'
                                        */
  real_T uL_Gain;                      /* Expression: 1/L
                                        * Referenced by: '<S85>/1//L'
                                        */
  real_T donotdeletethisgain_Gain_j;   /* Expression: 1
                                        * Referenced by: '<S19>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_i;   /* Expression: 1
                                        * Referenced by: '<S20>/do not delete this gain'
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
 * Block '<S5>/Add' : Unused code path elimination
 * Block '<S29>/0 1' : Unused code path elimination
 * Block '<S29>/Gain' : Unused code path elimination
 * Block '<S29>/Saturation' : Unused code path elimination
 * Block '<S29>/Sum' : Unused code path elimination
 * Block '<S29>/Switch' : Unused code path elimination
 * Block '<S29>/eee' : Unused code path elimination
 * Block '<S31>/0 1' : Unused code path elimination
 * Block '<S31>/1//Rsw' : Unused code path elimination
 * Block '<S31>/Switch' : Unused code path elimination
 * Block '<S6>/Add' : Unused code path elimination
 * Block '<S36>/0 1' : Unused code path elimination
 * Block '<S36>/Gain' : Unused code path elimination
 * Block '<S36>/Saturation' : Unused code path elimination
 * Block '<S36>/Sum' : Unused code path elimination
 * Block '<S36>/Switch' : Unused code path elimination
 * Block '<S36>/eee' : Unused code path elimination
 * Block '<S38>/0 1' : Unused code path elimination
 * Block '<S38>/1//Rsw' : Unused code path elimination
 * Block '<S38>/Switch' : Unused code path elimination
 * Block '<S7>/Add' : Unused code path elimination
 * Block '<S43>/0 1' : Unused code path elimination
 * Block '<S43>/Gain' : Unused code path elimination
 * Block '<S43>/Saturation' : Unused code path elimination
 * Block '<S43>/Sum' : Unused code path elimination
 * Block '<S43>/Switch' : Unused code path elimination
 * Block '<S43>/eee' : Unused code path elimination
 * Block '<S45>/0 1' : Unused code path elimination
 * Block '<S45>/1//Rsw' : Unused code path elimination
 * Block '<S45>/Switch' : Unused code path elimination
 * Block '<S8>/Add' : Unused code path elimination
 * Block '<S50>/0 1' : Unused code path elimination
 * Block '<S50>/Gain' : Unused code path elimination
 * Block '<S50>/Saturation' : Unused code path elimination
 * Block '<S50>/Sum' : Unused code path elimination
 * Block '<S50>/Switch' : Unused code path elimination
 * Block '<S50>/eee' : Unused code path elimination
 * Block '<S52>/0 1' : Unused code path elimination
 * Block '<S52>/1//Rsw' : Unused code path elimination
 * Block '<S52>/Switch' : Unused code path elimination
 * Block '<S9>/Add' : Unused code path elimination
 * Block '<S57>/0 1' : Unused code path elimination
 * Block '<S57>/Gain' : Unused code path elimination
 * Block '<S57>/Saturation' : Unused code path elimination
 * Block '<S57>/Sum' : Unused code path elimination
 * Block '<S57>/Switch' : Unused code path elimination
 * Block '<S57>/eee' : Unused code path elimination
 * Block '<S59>/0 1' : Unused code path elimination
 * Block '<S59>/1//Rsw' : Unused code path elimination
 * Block '<S59>/Switch' : Unused code path elimination
 * Block '<S10>/Add' : Unused code path elimination
 * Block '<S64>/0 1' : Unused code path elimination
 * Block '<S64>/Gain' : Unused code path elimination
 * Block '<S64>/Saturation' : Unused code path elimination
 * Block '<S64>/Sum' : Unused code path elimination
 * Block '<S64>/Switch' : Unused code path elimination
 * Block '<S64>/eee' : Unused code path elimination
 * Block '<S66>/0 1' : Unused code path elimination
 * Block '<S66>/1//Rsw' : Unused code path elimination
 * Block '<S66>/Switch' : Unused code path elimination
 * Block '<S11>/Add' : Unused code path elimination
 * Block '<S71>/0 1' : Unused code path elimination
 * Block '<S71>/Gain' : Unused code path elimination
 * Block '<S71>/Saturation' : Unused code path elimination
 * Block '<S71>/Sum' : Unused code path elimination
 * Block '<S71>/Switch' : Unused code path elimination
 * Block '<S71>/eee' : Unused code path elimination
 * Block '<S73>/0 1' : Unused code path elimination
 * Block '<S73>/1//Rsw' : Unused code path elimination
 * Block '<S73>/Switch' : Unused code path elimination
 * Block '<S12>/Add' : Unused code path elimination
 * Block '<S78>/0 1' : Unused code path elimination
 * Block '<S78>/Gain' : Unused code path elimination
 * Block '<S78>/Saturation' : Unused code path elimination
 * Block '<S78>/Sum' : Unused code path elimination
 * Block '<S78>/Switch' : Unused code path elimination
 * Block '<S78>/eee' : Unused code path elimination
 * Block '<S80>/0 1' : Unused code path elimination
 * Block '<S80>/1//Rsw' : Unused code path elimination
 * Block '<S80>/Switch' : Unused code path elimination
 * Block '<S13>/0' : Unused code path elimination
 * Block '<S14>/0' : Unused code path elimination
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
 * '<S5>'   : 'Simulacion_stepper/Mosfet'
 * '<S6>'   : 'Simulacion_stepper/Mosfet1'
 * '<S7>'   : 'Simulacion_stepper/Mosfet2'
 * '<S8>'   : 'Simulacion_stepper/Mosfet3'
 * '<S9>'   : 'Simulacion_stepper/Mosfet4'
 * '<S10>'  : 'Simulacion_stepper/Mosfet5'
 * '<S11>'  : 'Simulacion_stepper/Mosfet6'
 * '<S12>'  : 'Simulacion_stepper/Mosfet7'
 * '<S13>'  : 'Simulacion_stepper/Park to Clarke Angle Transform'
 * '<S14>'  : 'Simulacion_stepper/Park to Clarke Angle Transform1'
 * '<S15>'  : 'Simulacion_stepper/Stepper Motor'
 * '<S16>'  : 'Simulacion_stepper/Triangle Generator'
 * '<S17>'  : 'Simulacion_stepper/Triggered Subsystem'
 * '<S18>'  : 'Simulacion_stepper/Triggered Subsystem1'
 * '<S19>'  : 'Simulacion_stepper/Voltage Measurement'
 * '<S20>'  : 'Simulacion_stepper/Voltage Measurement1'
 * '<S21>'  : 'Simulacion_stepper/powergui'
 * '<S22>'  : 'Simulacion_stepper/Current Measurement/Model'
 * '<S23>'  : 'Simulacion_stepper/Current Measurement1/Model'
 * '<S24>'  : 'Simulacion_stepper/DC Voltage Source/Model'
 * '<S25>'  : 'Simulacion_stepper/DC Voltage Source1/Model'
 * '<S26>'  : 'Simulacion_stepper/Mosfet/Diode'
 * '<S27>'  : 'Simulacion_stepper/Mosfet/Ideal Switch'
 * '<S28>'  : 'Simulacion_stepper/Mosfet/Measurement list'
 * '<S29>'  : 'Simulacion_stepper/Mosfet/Diode/Model'
 * '<S30>'  : 'Simulacion_stepper/Mosfet/Diode/Model/Measurement list'
 * '<S31>'  : 'Simulacion_stepper/Mosfet/Ideal Switch/Model'
 * '<S32>'  : 'Simulacion_stepper/Mosfet/Ideal Switch/Model/Measurement list'
 * '<S33>'  : 'Simulacion_stepper/Mosfet1/Diode'
 * '<S34>'  : 'Simulacion_stepper/Mosfet1/Ideal Switch'
 * '<S35>'  : 'Simulacion_stepper/Mosfet1/Measurement list'
 * '<S36>'  : 'Simulacion_stepper/Mosfet1/Diode/Model'
 * '<S37>'  : 'Simulacion_stepper/Mosfet1/Diode/Model/Measurement list'
 * '<S38>'  : 'Simulacion_stepper/Mosfet1/Ideal Switch/Model'
 * '<S39>'  : 'Simulacion_stepper/Mosfet1/Ideal Switch/Model/Measurement list'
 * '<S40>'  : 'Simulacion_stepper/Mosfet2/Diode'
 * '<S41>'  : 'Simulacion_stepper/Mosfet2/Ideal Switch'
 * '<S42>'  : 'Simulacion_stepper/Mosfet2/Measurement list'
 * '<S43>'  : 'Simulacion_stepper/Mosfet2/Diode/Model'
 * '<S44>'  : 'Simulacion_stepper/Mosfet2/Diode/Model/Measurement list'
 * '<S45>'  : 'Simulacion_stepper/Mosfet2/Ideal Switch/Model'
 * '<S46>'  : 'Simulacion_stepper/Mosfet2/Ideal Switch/Model/Measurement list'
 * '<S47>'  : 'Simulacion_stepper/Mosfet3/Diode'
 * '<S48>'  : 'Simulacion_stepper/Mosfet3/Ideal Switch'
 * '<S49>'  : 'Simulacion_stepper/Mosfet3/Measurement list'
 * '<S50>'  : 'Simulacion_stepper/Mosfet3/Diode/Model'
 * '<S51>'  : 'Simulacion_stepper/Mosfet3/Diode/Model/Measurement list'
 * '<S52>'  : 'Simulacion_stepper/Mosfet3/Ideal Switch/Model'
 * '<S53>'  : 'Simulacion_stepper/Mosfet3/Ideal Switch/Model/Measurement list'
 * '<S54>'  : 'Simulacion_stepper/Mosfet4/Diode'
 * '<S55>'  : 'Simulacion_stepper/Mosfet4/Ideal Switch'
 * '<S56>'  : 'Simulacion_stepper/Mosfet4/Measurement list'
 * '<S57>'  : 'Simulacion_stepper/Mosfet4/Diode/Model'
 * '<S58>'  : 'Simulacion_stepper/Mosfet4/Diode/Model/Measurement list'
 * '<S59>'  : 'Simulacion_stepper/Mosfet4/Ideal Switch/Model'
 * '<S60>'  : 'Simulacion_stepper/Mosfet4/Ideal Switch/Model/Measurement list'
 * '<S61>'  : 'Simulacion_stepper/Mosfet5/Diode'
 * '<S62>'  : 'Simulacion_stepper/Mosfet5/Ideal Switch'
 * '<S63>'  : 'Simulacion_stepper/Mosfet5/Measurement list'
 * '<S64>'  : 'Simulacion_stepper/Mosfet5/Diode/Model'
 * '<S65>'  : 'Simulacion_stepper/Mosfet5/Diode/Model/Measurement list'
 * '<S66>'  : 'Simulacion_stepper/Mosfet5/Ideal Switch/Model'
 * '<S67>'  : 'Simulacion_stepper/Mosfet5/Ideal Switch/Model/Measurement list'
 * '<S68>'  : 'Simulacion_stepper/Mosfet6/Diode'
 * '<S69>'  : 'Simulacion_stepper/Mosfet6/Ideal Switch'
 * '<S70>'  : 'Simulacion_stepper/Mosfet6/Measurement list'
 * '<S71>'  : 'Simulacion_stepper/Mosfet6/Diode/Model'
 * '<S72>'  : 'Simulacion_stepper/Mosfet6/Diode/Model/Measurement list'
 * '<S73>'  : 'Simulacion_stepper/Mosfet6/Ideal Switch/Model'
 * '<S74>'  : 'Simulacion_stepper/Mosfet6/Ideal Switch/Model/Measurement list'
 * '<S75>'  : 'Simulacion_stepper/Mosfet7/Diode'
 * '<S76>'  : 'Simulacion_stepper/Mosfet7/Ideal Switch'
 * '<S77>'  : 'Simulacion_stepper/Mosfet7/Measurement list'
 * '<S78>'  : 'Simulacion_stepper/Mosfet7/Diode/Model'
 * '<S79>'  : 'Simulacion_stepper/Mosfet7/Diode/Model/Measurement list'
 * '<S80>'  : 'Simulacion_stepper/Mosfet7/Ideal Switch/Model'
 * '<S81>'  : 'Simulacion_stepper/Mosfet7/Ideal Switch/Model/Measurement list'
 * '<S82>'  : 'Simulacion_stepper/Stepper Motor/Model'
 * '<S83>'  : 'Simulacion_stepper/Stepper Motor/Model/FEM'
 * '<S84>'  : 'Simulacion_stepper/Stepper Motor/Model/Mechanical'
 * '<S85>'  : 'Simulacion_stepper/Stepper Motor/Model/windings'
 * '<S86>'  : 'Simulacion_stepper/Triangle Generator/Model'
 * '<S87>'  : 'Simulacion_stepper/Voltage Measurement/Model'
 * '<S88>'  : 'Simulacion_stepper/Voltage Measurement1/Model'
 * '<S89>'  : 'Simulacion_stepper/powergui/EquivalentModel1'
 * '<S90>'  : 'Simulacion_stepper/powergui/EquivalentModel1/Gates'
 * '<S91>'  : 'Simulacion_stepper/powergui/EquivalentModel1/Sources'
 * '<S92>'  : 'Simulacion_stepper/powergui/EquivalentModel1/Status'
 * '<S93>'  : 'Simulacion_stepper/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* Simulacion_stepper_h_ */
