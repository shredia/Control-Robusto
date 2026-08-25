function [Iq_comp, is_steady, LUT_mean, Wm_ripple, ...
          accel_out, t_steady_out] = ...
    online_lut_learning(Wm_filt, Wm_obs, theta_m, LUT_parameters)

% =========================================================================
% 0. PARÁMETROS
% =========================================================================

N_LUT         = LUT_parameters.N_LUT;
Wm_min        = LUT_parameters.Wm_min;
accel_tol     = LUT_parameters.accel_tol;
T_steady_req  = LUT_parameters.T_steady_req;

alpha_learn   = LUT_parameters.alpha_learn;
K_learn       = LUT_parameters.K_learn;

Ts            = LUT_parameters.Ts;
enable        = LUT_parameters.enable;


% =========================================================================
% PROTECCIONES
% =========================================================================

Ts = max(Ts, eps);

alpha_learn = ...
    min(max(alpha_learn, 0.0), 1.0);


% =========================================================================
% ESTADOS PERSISTENTES
% =========================================================================

persistent t_steady
persistent Wm_filt_prev
persistent accel_state
persistent LUT


if isempty(t_steady)

    t_steady     = 0.0;

    Wm_filt_prev = Wm_filt;

    accel_state  = 0.0;

    LUT = zeros(N_LUT,1);

end


% =========================================================================
% INICIALIZACIÓN DE SALIDAS
% =========================================================================

Iq_comp     = 0.0;
is_steady   = 0.0;
LUT_mean    = 0.0;
Wm_ripple   = 0.0;

accel_out    = 0.0;
t_steady_out = 0.0;


% =========================================================================
% 1. ACELERACIÓN DE Wm FILTRADA
% =========================================================================

accel_raw = ...
    (Wm_filt - Wm_filt_prev) / Ts;

Wm_filt_prev = Wm_filt;


% =========================================================================
% FILTRO DE ACELERACIÓN
% =========================================================================

tau_accel = 0.01;     % 10 ms

alpha_accel = ...
    Ts / (tau_accel + Ts);


accel_state = ...
    accel_state ...
    + alpha_accel * ...
      (accel_raw - accel_state);


accel_filt = accel_state;


% Salida para monitoreo
accel_out = accel_filt;


% =========================================================================
% 2. DETECCIÓN DE STEADY STATE
% =========================================================================

steady_condition = ...
       (abs(Wm_filt) >= Wm_min) ...
    && (abs(accel_filt) <= accel_tol);


if steady_condition

    t_steady = ...
        t_steady + Ts;

else

    t_steady = 0.0;

end


% Saturar contador
if t_steady >= T_steady_req

    t_steady = T_steady_req;

    is_steady = 1.0;

else

    is_steady = 0.0;

end


% Salida para monitoreo
t_steady_out = t_steady;


% =========================================================================
% 3. RIPPLE DE VELOCIDAD
% =========================================================================

Wm_ripple = ...
    Wm_obs - Wm_filt;


% =========================================================================
% 4. ÍNDICE ANGULAR LUT
% =========================================================================

theta_mod = ...
    mod(theta_m, 2*pi);


bin_idx = ...
    floor( ...
        theta_mod / (2*pi) * N_LUT ...
    ) + 1;


if bin_idx > N_LUT

    bin_idx = N_LUT;

elseif bin_idx < 1

    bin_idx = 1;

end


% =========================================================================
% 5. APRENDIZAJE
% =========================================================================

if (enable > 0.5) && ...
   (is_steady > 0.5)

    delta_Iq = ...
        -K_learn * Wm_ripple;


    LUT(bin_idx) = ...
        LUT(bin_idx) ...
        + alpha_learn * delta_Iq;

end


% =========================================================================
% 6. COMPENSACIÓN
% =========================================================================

if enable > 0.5

    Iq_comp = LUT(bin_idx);

else

    Iq_comp = 0.0;

end


% =========================================================================
% 7. MONITOREO LUT
% =========================================================================

LUT_mean = ...
    sum(LUT) / N_LUT;


end