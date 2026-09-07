%% ========================================================================
% INITIALIZATION - CONTROL ROBUSTO
% ========================================================================

project_folder = fileparts(mfilename('fullpath'));

% Carpeta principal
addpath(project_folder);

% Controladores
addpath(fullfile(project_folder,'control'));

fprintf('Proyecto inicializado desde:\n%s\n', project_folder);


%% ========================================================================
% SIMULACIÓN Y REFERENCIAS
% ========================================================================

Time_simulation = 0.2;

t_ref = [0.0];
w_ref = [10];

theta_ref = pi/4;

Wm_ref = timeseries(w_ref, t_ref);
Wm_ref = setinterpmethod(Wm_ref,'zoh');


% -------------------------------------------------------------------------
% Torque de carga
% -------------------------------------------------------------------------

t_load = [0.0 0.5];
Tl_ref = [0.0 0.0];

Tl = timeseries(Tl_ref, t_load);
Tl = setinterpmethod(Tl,'linear');


% -------------------------------------------------------------------------
% Tensión DC
% -------------------------------------------------------------------------

Vdc = 24;


%% ========================================================================
% TIEMPOS DE MUESTREO Y FRECUENCIAS
% ========================================================================

frecuency_simulation = 10e3;

sample_time = ...
    1/frecuency_simulation;


% -------------------------------------------------------------------------
% PWM
% -------------------------------------------------------------------------

f_carrier = 20e3;


% -------------------------------------------------------------------------
% Control de corriente
% -------------------------------------------------------------------------

f_current = ...
    f_carrier/2;

Ts_current = ...
    1/f_current;


% -------------------------------------------------------------------------
% EKF
% -------------------------------------------------------------------------

f_ekf = ...
    f_current;

Ts_ekf = ...
    1/f_ekf;


% -------------------------------------------------------------------------
% Control de velocidad
% -------------------------------------------------------------------------

Ts_Wm = ...
    Ts_current*10;

f_wm = ...
    1/Ts_Wm;


% -------------------------------------------------------------------------
% Control de posición
% -------------------------------------------------------------------------

Ts_pos = ...
    Ts_Wm*10;


% -------------------------------------------------------------------------
% Disturbance Observer
% -------------------------------------------------------------------------

Ts_DO = ...
    Ts_Wm;


%% ========================================================================
% PARÁMETROS REALES DEL MOTOR
% ========================================================================

motor_params.Step_angle = 1.8;

motor_params.N_phases = 2;


% -------------------------------------------------------------------------
% Geometría
% -------------------------------------------------------------------------

motor_params.N_steps = ...
    360/motor_params.Step_angle;


motor_params.N_teeths = ...
    motor_params.N_steps / ...
    motor_params.N_phases;


motor_params.P = ...
    motor_params.N_teeths/2;


motor_params.Nr = ...
    motor_params.N_steps;


% -------------------------------------------------------------------------
% Corrientes
% -------------------------------------------------------------------------

motor_params.I_nom = 1.7;


motor_params.Imax = ...
    sqrt(2)*motor_params.I_nom;


% -------------------------------------------------------------------------
% Torque
% -------------------------------------------------------------------------

motor_params.Tdm = ...
    22/1000;


motor_params.Thold = ...
    392/1000;


% -------------------------------------------------------------------------
% Constantes electromecánicas
% -------------------------------------------------------------------------

motor_params.Kt = ...
    motor_params.Thold / ...
    (sqrt(2)*motor_params.I_nom);


motor_params.Psi = ...
    motor_params.Kt / ...
    motor_params.P;


motor_params.Ke = ...
    motor_params.Kt;


% -------------------------------------------------------------------------
% Parámetros eléctricos
% -------------------------------------------------------------------------

motor_params.R = 2.5;

motor_params.L = 6/1000;

motor_params.Lq = 9.61/1000;

motor_params.Ld = 3.66/1000;

motor_params.L0 = 6.635/1000;

motor_params.L2 = 2.975/1000;


% -------------------------------------------------------------------------
% Fase del cogging real
% -------------------------------------------------------------------------

motor_params.Phi = ...
    pi/2;


% -------------------------------------------------------------------------
% Parámetros mecánicos internos
% -------------------------------------------------------------------------

motor_params.B_internal = ...
    1e-4;


motor_params.J_internal = ...
    54e-7;


% -------------------------------------------------------------------------
% Carga externa
% -------------------------------------------------------------------------

motor_params.J_external = ...
    motor_params.J_internal*0;


motor_params.B_external = ...
    motor_params.B_internal*0;


% -------------------------------------------------------------------------
% Parámetros mecánicos reales totales
% -------------------------------------------------------------------------

motor_params.J_real = ...
    motor_params.J_internal + ...
    motor_params.J_external;


motor_params.B_real = ...
    motor_params.B_internal + ...
    motor_params.B_external;


% -------------------------------------------------------------------------
% Fricción
% -------------------------------------------------------------------------

motor_params.Tc = 0.002;

motor_params.Tba = 0;

motor_params.Wba = 0;


% -------------------------------------------------------------------------
% Condiciones iniciales
% -------------------------------------------------------------------------

motor_params.InitialSpeed = 0;


% -------------------------------------------------------------------------
% Otros
% -------------------------------------------------------------------------

motor_params.max_step_rate = 3000;


%% ========================================================================
% CONTROL DE VELOCIDAD + MTPA
% ========================================================================

PI_MTPA_params.fbw_Wm = 5;

PI_MTPA_params.shi_w = 2;


PI_MTPA_params.wn_w = ...
    2*pi*PI_MTPA_params.fbw_Wm;


% -------------------------------------------------------------------------
% Parámetros mecánicos conocidos por el controlador
% -------------------------------------------------------------------------

PI_MTPA_params.J = ...
    motor_params.J_internal;


PI_MTPA_params.B = ...
    motor_params.B_internal;


% -------------------------------------------------------------------------
% Ganancias PI/PID de velocidad
% -------------------------------------------------------------------------

PI_MTPA_params.Kp_w = ...
    ( ...
        2 * ...
        PI_MTPA_params.shi_w * ...
        PI_MTPA_params.wn_w * ...
        PI_MTPA_params.J ...
    ) / ...
    motor_params.Kt;


PI_MTPA_params.Ki_w = ...
    (PI_MTPA_params.wn_w^2) * ...
    PI_MTPA_params.J / ...
    motor_params.Kt;


PI_MTPA_params.Kd_w = ...
    PI_MTPA_params.Kp_w/10000;


PI_MTPA_params.Tf_w = ...
    10/1000;


% -------------------------------------------------------------------------
% Feedforward
% -------------------------------------------------------------------------

PI_MTPA_params.Tl_tdm = ...
    motor_params.Tdm;


% -------------------------------------------------------------------------
% Límite de corriente
% -------------------------------------------------------------------------

PI_MTPA_params.Imax = ...
    motor_params.Imax;


% -------------------------------------------------------------------------
% Parámetros eléctricos usados por MTPA
% -------------------------------------------------------------------------

PI_MTPA_params.Ld = ...
    motor_params.Ld;


PI_MTPA_params.Lq = ...
    motor_params.Lq;


PI_MTPA_params.Ke = ...
    motor_params.Ke;


PI_MTPA_params.P = ...
    motor_params.P;


PI_MTPA_params.Tmax = ...
    0.5;


% -------------------------------------------------------------------------
% Tiempo de muestreo
% -------------------------------------------------------------------------

PI_MTPA_params.Ts = ...
    Ts_Wm;


%% ========================================================================
% CONTROL DE POSICIÓN
% ========================================================================

PI_pos_params.Kp = 0.5;

PI_pos_params.Ki = 0.2;

PI_pos_params.Kd = 0.001;

PI_pos_params.Ts = ...
    Ts_pos;


PI_position_params.Ts = ...
    Ts_current;


PI_position_params.Wmax_pos = ...
    20;


PI_position_params.Amax_pos = ...
    100;


PI_position_params.theta_tol = ...
    0.01;


%% ========================================================================
% PARÁMETROS DEL EKF
% ========================================================================

EKF_params.Ts = ...
    Ts_ekf;


% -------------------------------------------------------------------------
% Modelo eléctrico
% -------------------------------------------------------------------------

EKF_params.R = ...
    motor_params.R;


EKF_params.Ld = ...
    motor_params.Ld;


EKF_params.Lq = ...
    motor_params.Lq;


EKF_params.Ke = ...
    motor_params.Ke;


EKF_params.Kt = ...
    motor_params.Kt;


EKF_params.P = ...
    motor_params.P;


% -------------------------------------------------------------------------
% Modelo mecánico
% -------------------------------------------------------------------------

EKF_params.J = ...
    motor_params.J_internal;


EKF_params.B = ...
    motor_params.B_internal;


% -------------------------------------------------------------------------
% Cogging
% -------------------------------------------------------------------------

EKF_params.Tdm = ...
    motor_params.Tdm;


EKF_params.N_steps = ...
    motor_params.N_steps;


EKF_params.Nr = ...
    motor_params.Nr;


EKF_params.Phi = ...
    motor_params.Phi;


% -------------------------------------------------------------------------
% Límites de estados
% -------------------------------------------------------------------------

EKF_params.Id_max = 3;

EKF_params.Iq_max = 3;

EKF_params.Wm_max = 40;

EKF_params.Tx_max = 0.5;


%% ========================================================================
% HFI
% ========================================================================

EKF_params.HFI_flag = 0;


EKF_params.Amplitud_HFI = 2;


EKF_params.f_h = 2000;


EKF_params.wh = ...
    2*pi*EKF_params.f_h;


EKF_params.BW_hfi = 25;


Fs_control = ...
    f_ekf;


Wn = [ ...
    (EKF_params.f_h - EKF_params.BW_hfi/2)/(Fs_control/2), ...
    (EKF_params.f_h + EKF_params.BW_hfi/2)/(Fs_control/2) ...
];


[b_hfi,a_hfi] = ...
    butter( ...
        1, ...
        Wn, ...
        'bandpass');


EKF_params.b0_hfi = ...
    b_hfi(1);


EKF_params.b1_hfi = ...
    b_hfi(2);


EKF_params.b2_hfi = ...
    b_hfi(3);


EKF_params.a1_hfi = ...
    a_hfi(2);


EKF_params.a2_hfi = ...
    a_hfi(3);


fprintf( ...
    'Coeficientes HFI para f_h = %d Hz:\n', ...
    EKF_params.f_h);


fprintf( ...
    'b0 = %.6f, b1 = %.6f, b2 = %.6f\n', ...
    EKF_params.b0_hfi, ...
    EKF_params.b1_hfi, ...
    EKF_params.b2_hfi);


fprintf( ...
    'a1 = %.6f, a2 = %.6f\n', ...
    EKF_params.a1_hfi, ...
    EKF_params.a2_hfi);


%% ========================================================================
% CONTROLADOR DE CORRIENTE PI dq
% ========================================================================

PI_dq_params.Vdc = ...
    Vdc;


% -------------------------------------------------------------------------
% Ancho de banda
% -------------------------------------------------------------------------

PI_dq_params.fbw_d = 500;

PI_dq_params.fbw_q = 500;


PI_dq_params.wd_d = ...
    2*pi*PI_dq_params.fbw_d;


PI_dq_params.wd_q = ...
    2*pi*PI_dq_params.fbw_q;


PI_dq_params.shi_d = 1;

PI_dq_params.shi_q = 1;


% -------------------------------------------------------------------------
% Control usando inductancia nominal
% -------------------------------------------------------------------------

PI_dq_params.Kp_d = ...
    PI_dq_params.wd_d * ...
    motor_params.L;


PI_dq_params.Ki_d = ...
    PI_dq_params.wd_d * ...
    motor_params.R;


PI_dq_params.Kp_q = ...
    PI_dq_params.wd_q * ...
    motor_params.L;


PI_dq_params.Ki_q = ...
    PI_dq_params.wd_q * ...
    motor_params.R;


% -------------------------------------------------------------------------
% Control considerando saliencia
% -------------------------------------------------------------------------

PI_dq_params.Kp_q_salient = ...
    2 * ...
    PI_dq_params.shi_q * ...
    PI_dq_params.wd_q * ...
    motor_params.Lq ...
    - motor_params.R;


PI_dq_params.Ki_q_salient = ...
    PI_dq_params.wd_q^2 * ...
    motor_params.Lq;


PI_dq_params.Kp_d_salient = ...
    2 * ...
    PI_dq_params.shi_d * ...
    PI_dq_params.wd_d * ...
    motor_params.Ld ...
    - motor_params.R;


PI_dq_params.Ki_d_salient = ...
    PI_dq_params.wd_d^2 * ...
    motor_params.Ld;


% -------------------------------------------------------------------------
% Tiempo de muestreo
% -------------------------------------------------------------------------

PI_dq_params.Ts = ...
    Ts_current;


% -------------------------------------------------------------------------
% Modelo
% -------------------------------------------------------------------------

PI_dq_params.Ld = ...
    motor_params.Ld;


PI_dq_params.Lq = ...
    motor_params.Lq;


PI_dq_params.Psi = ...
    motor_params.Psi;


PI_dq_params.P = ...
    motor_params.P;


PI_dq_params.Imax = ...
    motor_params.Imax;


% -------------------------------------------------------------------------
% HFI
% -------------------------------------------------------------------------

PI_dq_params.Amplitud_HFI = ...
    EKF_params.Amplitud_HFI;


PI_dq_params.wh = ...
    EKF_params.wh;


PI_dq_params.HFI_flag = ...
    EKF_params.HFI_flag;


%% ========================================================================
% CONTROL RESONANTE DE VELOCIDAD
%
% M1:
%
%   Wm_obs ---> resonador ---> T_PR_M1
%
%   Utiliza:
%
%       RI_params_M1
%
%   La fase puede ser compensada:
%
%       RI_params_M1.phi
%
%
% M2:
%
%   Wm_real ---> resonador ---> T_PR_M2
%
%   Utiliza:
%
%       RI_params_M2
%
%   Se deja inicialmente:
%
%       phi_M2 = 0 deg
%
% =========================================================================


%% ========================================================================
% PARÁMETROS BASE
% ========================================================================

RI_base.enable = 1;


% -------------------------------------------------------------------------
% Tiempo de muestreo
% -------------------------------------------------------------------------

RI_base.Ts = ...
    Ts_current;


% -------------------------------------------------------------------------
% Número de ciclos de cogging por revolución
% -------------------------------------------------------------------------

RI_base.Nc = ...
    200;


% -------------------------------------------------------------------------
% Ganancia resonante
% -------------------------------------------------------------------------

RI_base.Kr = ...
    0.05;


% -------------------------------------------------------------------------
% Ancho del resonador
% -------------------------------------------------------------------------

RI_base.wc = ...
    5;


% -------------------------------------------------------------------------
% Saturación del torque resonante
% -------------------------------------------------------------------------

RI_base.T_res_max = ...
    0.021;


% -------------------------------------------------------------------------
% Límites de frecuencia
% -------------------------------------------------------------------------

RI_base.wr_min = ...
    10;


RI_base.wr_max = ...
    0.8*pi/Ts_current;


%% ========================================================================
% FASE DEL MOTOR 1
%
% M1 = caso sensorless
% =========================================================================

if ~exist('phi_M1_deg','var')

    phi_M1_deg = 0;

end


%% ========================================================================
% FASE DEL MOTOR 2
%
% M2 = referencia utilizando velocidad real
%
% Por defecto:
%
%       phi_M2 = 0 deg
%
% =========================================================================

if ~exist('phi_M2_deg','var')

    phi_M2_deg = 0;

end


%% ========================================================================
% CREAR PARÁMETROS INDEPENDIENTES
% ========================================================================

RI_params_M1 = ...
    RI_base;


RI_params_M2 = ...
    RI_base;


%% ========================================================================
% ASIGNAR FASES
% ========================================================================

RI_params_M1.phi = ...
    phi_M1_deg*pi/180;


RI_params_M2.phi = ...
    phi_M2_deg*pi/180;


%% ========================================================================
% INFORMACIÓN DEL CONTROL RESONANTE
% ========================================================================

fprintf('\n');

fprintf('=============================================================\n');
fprintf(' CONTROL RESONANTE DE VELOCIDAD\n');
fprintf('=============================================================\n');


fprintf( ...
    'M1 - Wm observada : phi = %+8.3f deg\n', ...
    phi_M1_deg);


fprintf( ...
    'M2 - Wm real      : phi = %+8.3f deg\n', ...
    phi_M2_deg);


fprintf( ...
    'Kr               : %.6f\n', ...
    RI_base.Kr);


fprintf( ...
    'wc               : %.6f rad/s\n', ...
    RI_base.wc);


fprintf( ...
    'Nc               : %d\n', ...
    RI_base.Nc);


fprintf( ...
    'T_res max        : %.6f Nm\n', ...
    RI_base.T_res_max);


fprintf('=============================================================\n');


%% ========================================================================
% PR DE CORRIENTE d-q
%
% Este resonador es INDEPENDIENTE del resonador de velocidad M1/M2.
%
% La frecuencia sigue:
%
%       wr = Nc * abs(Wm)
%
% =========================================================================


% -------------------------------------------------------------------------
% Activación
%
% IMPORTANTE:
%
% Antes:
%
%       PR_i.enable = RI_params.enable;
%
% Ahora RI_params ya no existe.
%
% -------------------------------------------------------------------------

PR_i.enable = ...
    RI_base.enable;


% -------------------------------------------------------------------------
% Número de ciclos de cogging por revolución
% -------------------------------------------------------------------------

PR_i.Nc = ...
    200;


% -------------------------------------------------------------------------
% Ancho del resonador
% -------------------------------------------------------------------------

PR_i.wc = ...
    5;


% -------------------------------------------------------------------------
% Ganancias resonantes
% -------------------------------------------------------------------------

PR_i.Kr_d = ...
    0;


PR_i.Kr_q = ...
    1;


% -------------------------------------------------------------------------
% Saturación de los resonadores de corriente [V]
% -------------------------------------------------------------------------

PR_i.Vres_max_d = ...
    1.0;


PR_i.Vres_max_q = ...
    1.0;


% -------------------------------------------------------------------------
% Frecuencia resonante mínima
% -------------------------------------------------------------------------

PR_i.wr_min = ...
    10;


% -------------------------------------------------------------------------
% Frecuencia resonante máxima
% -------------------------------------------------------------------------

PR_i.wr_max = ...
    0.8*pi/Ts_current;


% -------------------------------------------------------------------------
% Guardar dentro de PI_dq
% -------------------------------------------------------------------------

PI_dq_params.PR_i = ...
    PR_i;


%% ========================================================================
% RESUMEN FINAL
% ========================================================================

fprintf('\n');

fprintf('=============================================================\n');
fprintf(' RESUMEN DE INICIALIZACIÓN\n');
fprintf('=============================================================\n');


fprintf( ...
    'Wm referencia       : %.4f rad/s\n', ...
    w_ref(end));


fprintf( ...
    'Ts corriente        : %.2e s\n', ...
    Ts_current);


fprintf( ...
    'Ts EKF              : %.2e s\n', ...
    Ts_ekf);


fprintf( ...
    'Ts velocidad        : %.2e s\n', ...
    Ts_Wm);


fprintf( ...
    'Cogging Nc          : %d\n', ...
    RI_base.Nc);


fprintf( ...
    'fcog @ Wm_ref       : %.4f Hz\n', ...
    RI_base.Nc*abs(w_ref(end))/(2*pi));


fprintf( ...
    'Phi M1              : %+8.3f deg\n', ...
    phi_M1_deg);


fprintf( ...
    'Phi M2              : %+8.3f deg\n', ...
    phi_M2_deg);


fprintf('=============================================================\n');