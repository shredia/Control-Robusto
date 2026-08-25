project_folder = fileparts(get_param(bdroot,'FileName'));

addpath(project_folder);
addpath(fullfile(project_folder,'control'));

%% =========================================================================
% SIMULACIÓN Y REFERENCIAS
% =========================================================================

Time_simulation = 0.5;

t_ref = [0.0];
w_ref = [10];

theta_ref = pi/4;

Wm_ref = timeseries(w_ref, t_ref);
Wm_ref = setinterpmethod(Wm_ref,'zoh');

% Torque de carga
t_load = [0.0 0.5];
Tl_ref = [0.0 0.0];

Tl = timeseries(Tl_ref, t_load);
Tl = setinterpmethod(Tl, 'linear');

Vdc = 24;

%% =========================================================================
% TIEMPOS DE MUESTREO Y FRECUENCIAS
% =========================================================================

frecuency_simulation = 10e3;
sample_time = 1/frecuency_simulation;

f_carrier = 20e3;

f_current  = f_carrier/2;
Ts_current = 1/f_current;

f_ekf  = f_current;
Ts_ekf = 1/f_ekf;

Ts_Wm  = Ts_current*10;
Ts_pos = Ts_Wm*10;

f_wm = 1/Ts_Wm;

Ts_DO = Ts_Wm;



%% =========================================================================
% PARÁMETROS REALES DEL MOTOR
% =========================================================================

motor_params.Step_angle = 1.8;
motor_params.N_phases   = 2;

motor_params.N_steps = ...
    360/motor_params.Step_angle;

motor_params.N_teeths = ...
    motor_params.N_steps/motor_params.N_phases;

motor_params.P = ...
    motor_params.N_teeths/2;

% Corrientes
motor_params.I_nom = 1.7;

motor_params.Imax = ...
    sqrt(2)*motor_params.I_nom;

% Torque
motor_params.Tdm   = 22/1000;
motor_params.Thold = 392/1000;

% Constantes electromecánicas
motor_params.Kt = ...
    motor_params.Thold / ...
    (sqrt(2)*motor_params.I_nom);

motor_params.Psi = ...
    motor_params.Kt/motor_params.P;

motor_params.Ke = motor_params.Kt;

% Parámetros eléctricos
motor_params.R = 2.5;

motor_params.L = 6/1000;

motor_params.Lq = 9.61/1000;
motor_params.Ld = 3.66/1000;

motor_params.L0 = 6.635/1000;
motor_params.L2 = 2.975/1000;

motor_params.Phi = pi/2;

% Mecánicos internos
motor_params.B_internal = 1e-4;
motor_params.J_internal = 54e-7;

% Carga externa
motor_params.J_external = ...
    motor_params.J_internal*0;

motor_params.B_external = ...
    motor_params.B_internal*0;

% Parámetros reales totales
motor_params.J_real = ...
    motor_params.J_internal + ...
    motor_params.J_external;

motor_params.B_real = ...
    motor_params.B_internal + ...
    motor_params.B_external;

% Fricción
motor_params.Tc  = 0.002;
motor_params.Tba = 0;
motor_params.Wba = 0;

motor_params.InitialSpeed = 0;

% Otros
motor_params.max_step_rate = 3000;

motor_params.Nr = ...
    motor_params.N_steps;




%% =========================================================================
% PARÁMETROS DEL CONTROL DE VELOCIDAD + MTPA
% =========================================================================

PI_MTPA_params.fbw_Wm = 10;%%10 para control de velocidad.

PI_MTPA_params.shi_w = 2;

PI_MTPA_params.wn_w = ...
    2*pi*PI_MTPA_params.fbw_Wm;

% Parámetros "conocidos" o asumidos por el controlador
PI_MTPA_params.J = motor_params.J_internal;
PI_MTPA_params.B = motor_params.B_internal;

PI_MTPA_params.Kp_w = ...
    (2 * PI_MTPA_params.shi_w * ...
     PI_MTPA_params.wn_w * ...
     PI_MTPA_params.J) / ...
     motor_params.Kt;

PI_MTPA_params.Ki_w = ...
    (PI_MTPA_params.wn_w^2) * ...
    PI_MTPA_params.J / ...
    motor_params.Kt;

PI_MTPA_params.Kd_w = ...
    PI_MTPA_params.Kp_w/10000;

PI_MTPA_params.Tf_w = 10/1000;

% Feedforward
PI_MTPA_params.Tl_tdm = ...
    motor_params.Tdm;

PI_MTPA_params.B = ...
    motor_params.B_internal;

% Límite de corriente
PI_MTPA_params.Imax = ...
    motor_params.Imax;

% Parámetros eléctricos usados por MTPA
PI_MTPA_params.Ld = motor_params.Ld;
PI_MTPA_params.Lq = motor_params.Lq;
PI_MTPA_params.Ke = motor_params.Ke;
PI_MTPA_params.P  = motor_params.P;
PI_MTPA_params.Tmax = 0.5;



PI_MTPA_params.Ts = Ts_Wm;


%% =========================================================================
% CONTROL DE POSICIÓN
% =========================================================================

PI_pos_params.Kp = 0.5;
PI_pos_params.Ki = 0.2;
PI_pos_params.Kd = 0.001;

PI_pos_params.Ts = Ts_pos;

PI_position_params.Ts = Ts_current;

PI_position_params.Wmax_pos = 20;       % rad/s
PI_position_params.Amax_pos = 100;      % rad/s^2

PI_position_params.theta_tol = 0.01;    % rad


%% =========================================================================
% PARÁMETROS DEL EKF
% =========================================================================

EKF_params.Ts = Ts_ekf;

% Parámetros del modelo usado por el estimador
EKF_params.R  = motor_params.R;
EKF_params.Ld = motor_params.Ld;
EKF_params.Lq = motor_params.Lq;

EKF_params.Ke = motor_params.Ke;
EKF_params.Kt = motor_params.Kt;

EKF_params.P = motor_params.P;

% Aquí puedes colocar valores distintos a los reales
EKF_params.J = motor_params.J_internal;
EKF_params.B = motor_params.B_internal;

EKF_params.Tdm = motor_params.Tdm;

EKF_params.N_steps = motor_params.N_steps;
EKF_params.Nr      = motor_params.Nr;
EKF_params.Id_max = 3;
EKF_params.Iq_max = 3;
EKF_params.Wm_max = 40;
EKF_params.Tx_max = 0.5;



%% =========================================================================
% HFI
% =========================================================================

EKF_params.HFI_flag = 1;

EKF_params.Amplitud_HFI = 2;

EKF_params.f_h = 1000;

EKF_params.wh = ...
    2*pi*EKF_params.f_h;

EKF_params.BW_hfi = 100;

Fs_control = f_ekf;

Wn = [ ...
    (EKF_params.f_h - EKF_params.BW_hfi/2)/(Fs_control/2), ...
    (EKF_params.f_h + EKF_params.BW_hfi/2)/(Fs_control/2)];

[b_hfi, a_hfi] = ...
    butter(1, Wn, 'bandpass');

EKF_params.b0_hfi = b_hfi(1);
EKF_params.b1_hfi = b_hfi(2);
EKF_params.b2_hfi = b_hfi(3);

EKF_params.a1_hfi = a_hfi(2);
EKF_params.a2_hfi = a_hfi(3);

fprintf('Coeficientes HFI para f_h = %d Hz:\n', ...
    EKF_params.f_h);

fprintf('b0 = %.6f, b1 = %.6f, b2 = %.6f\n', ...
    EKF_params.b0_hfi, ...
    EKF_params.b1_hfi, ...
    EKF_params.b2_hfi);

fprintf('a1 = %.6f, a2 = %.6f\n', ...
    EKF_params.a1_hfi, ...
    EKF_params.a2_hfi);

%% =========================================================================
% PARÁMETROS DEL CONTROLADOR DE CORRIENTE PI dq
% =========================================================================

PI_dq_params.Vdc = Vdc;

PI_dq_params.fbw_d = 500;
PI_dq_params.fbw_q = 500;

PI_dq_params.wd_d = ...
    2*pi*PI_dq_params.fbw_d;

PI_dq_params.wd_q = ...
    2*pi*PI_dq_params.fbw_q;

PI_dq_params.shi_d = 1;
PI_dq_params.shi_q = 1;

% Control usando inductancia nominal
PI_dq_params.Kp_d = ...
    PI_dq_params.wd_d * motor_params.L;

PI_dq_params.Ki_d = ...
    PI_dq_params.wd_d * motor_params.R;

PI_dq_params.Kp_q = ...
    PI_dq_params.wd_q * motor_params.L;

PI_dq_params.Ki_q = ...
    PI_dq_params.wd_q * motor_params.R;

% Control considerando saliencia
PI_dq_params.Kp_q_salient = ...
    2*PI_dq_params.shi_q* ...
    PI_dq_params.wd_q* ...
    motor_params.Lq ...
    - motor_params.R;

PI_dq_params.Ki_q_salient = ...
    PI_dq_params.wd_q^2 * ...
    motor_params.Lq;

PI_dq_params.Kp_d_salient = ...
    2*PI_dq_params.shi_d* ...
    PI_dq_params.wd_d* ...
    motor_params.Ld ...
    - motor_params.R;

PI_dq_params.Ki_d_salient = ...
    PI_dq_params.wd_d^2 * ...
    motor_params.Ld;

PI_dq_params.Ts = Ts_current;

PI_dq_params.Ld = motor_params.Ld;
PI_dq_params.Lq = motor_params.Lq;
PI_dq_params.Psi = motor_params.Psi;
PI_dq_params.P = motor_params.P;
PI_dq_params.Imax = motor_params.Imax;
PI_dq_params.Amplitud_HFI = EKF_params.Amplitud_HFI;
PI_dq_params.wh = EKF_params.wh;
PI_dq_params.HFI_flag = EKF_params.HFI_flag;



%% ========================================================================
% LUT ONLINE COGGING
% =========================================================================

LUT_parameters.N_LUT         = 180;

LUT_parameters.Wm_min        = 1.0;
LUT_parameters.accel_tol     = 10;
LUT_parameters.T_steady_req  = 0.05;

LUT_parameters.alpha_learn   = 0.1;
LUT_parameters.K_learn       = 0.5;

LUT_parameters.Iq_comp_max   = 0.20;

LUT_parameters.Ts            = Ts_current;
LUT_parameters.enable        = 1;



% =========================================================================
% CONTROL RESONANTE
% =========================================================================

RI_params.enable = 0;

RI_params.Ts = Ts_Wm;

% Velocidad fija de operación
Wm_operacion = Wm_ref.Data;        % [rad/s]

% Número de ciclos de cogging por revolución
Nc = 200;                  % verificar según Tdtm real

% Frecuencia resonante
RI_params.wr =  1000*2*pi;    % [rad/s]

% Ganancia del resonador
RI_params.Kr = 0.05;

% Amortiguamiento del resonador
RI_params.zeta = 0.05;

% Torque máximo permitido al resonador
RI_params.T_res_max = 0.07;    % [Nm]

PI_MTPA_params.RI = RI_params;
