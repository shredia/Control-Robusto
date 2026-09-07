function [T_ref, T_PI, error_Wm, Wm_filt] = ...
    PI_speed(Wm_ref, Wm, Tl_tdm, params)

% =========================================================================
% CONTROL DE VELOCIDAD PI/PID
%
% La velocidad medida Wm se filtra mediante un LPF antes de entrar
% al controlador.
%
% Arquitectura:
%
%                Wm
%                 |
%                 v
%                LPF
%                 |
%                 v
%              Wm_filt
%                 |
%                 v
%       error = Wm_ref - Wm_filt
%                 |
%                 v
%              PI/PID
%                 |
%                T_PI
%                 |
%              + T_ff
%                 |
%                 v
%               T_ref
%
% =========================================================================


% =========================================================================
% 0. Parámetros
% =========================================================================

Ts = params.Ts;

Kp_w = params.Kp_w;
Ki_w = params.Ki_w;
Kd_w = params.Kd_w;

Tf_w = params.Tf_w;

B = params.B;

Tmax = params.Tmax;


% =========================================================================
% 1. Estados persistentes
% =========================================================================

persistent int_w
persistent Wm_filt_prev
persistent Wm_filt_der_prev
persistent dWm_filt_prev


if isempty(int_w)

    int_w = 0.0;

    Wm_filt_prev = 0.0;

    Wm_filt_der_prev = 0.0;

    dWm_filt_prev = 0.0;

end


% =========================================================================
% 2. Protección numérica
% =========================================================================

Ts = max(Ts,1e-9);

Tf_w = max(Tf_w,Ts);

Tmax = max(abs(Tmax),1e-9);


% =========================================================================
% 3. FILTRO PASABAJOS DE VELOCIDAD
% =========================================================================

alpha_w = ...
    Ts/(Tf_w + Ts);


Wm_filt = ...
    Wm_filt_prev ...
    + alpha_w * ...
    (Wm - Wm_filt_prev);


% =========================================================================
% 4. ERROR DE VELOCIDAD
% =========================================================================

error_Wm = ...
    Wm_ref - Wm; %%Ignoramos el filtro


% =========================================================================
% 5. DERIVADA DE VELOCIDAD FILTRADA
% =========================================================================

dWm_raw = ...
    (Wm_filt - Wm_filt_der_prev)/Ts;


% =========================================================================
% 6. FILTRO DE LA DERIVADA
% =========================================================================

alpha_d = ...
    Ts/(Tf_w + Ts);


dWm_filt = ...
    dWm_filt_prev ...
    + alpha_d * ...
    (dWm_raw - dWm_filt_prev);


% =========================================================================
% 7. TÉRMINO PROPORCIONAL
% =========================================================================

T_P = ...
    Kp_w * error_Wm;


% =========================================================================
% 8. TÉRMINO DERIVATIVO
%
% Derivada sobre la medida para evitar derivative kick.
% =========================================================================

T_D = ...
    -Kd_w * dWm_filt;


% =========================================================================
% 9. TORQUE PI/PID
% =========================================================================

T_PI = ...
    T_P ...
    + int_w ...
    + T_D;


% =========================================================================
% 10. FEEDFORWARD
%
% Compensación del torque de carga estimado y fricción viscosa.
% =========================================================================

T_ff = ...
    Tl_tdm ...
    + B*Wm_filt;


% =========================================================================
% 11. TORQUE TOTAL SIN SATURAR
% =========================================================================

T_unsat = ...
    T_PI ...
    + T_ff;


% =========================================================================
% 12. SATURACIÓN DE TORQUE
% =========================================================================

T_ref = ...
    min( ...
        max(T_unsat,-Tmax), ...
        Tmax);


% =========================================================================
% 13. ANTI-WINDUP
%
% Integración condicional:
%
% - integra normalmente si no existe saturación
% - si está saturado, solo integra cuando el error ayuda a salir
%   de la saturación
% =========================================================================

sat_high = ...
    T_unsat > Tmax;

sat_low = ...
    T_unsat < -Tmax;


integrate = ...
       (~sat_high && ~sat_low) ...
    || (sat_high && error_Wm < 0.0) ...
    || (sat_low  && error_Wm > 0.0);


if integrate

    int_w = ...
        int_w ...
        + Ki_w*error_Wm*Ts;

end


% =========================================================================
% 14. ACTUALIZAR ESTADOS
% =========================================================================

Wm_filt_prev = Wm_filt;

Wm_filt_der_prev = Wm_filt;

dWm_filt_prev = dWm_filt;


end