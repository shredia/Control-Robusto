function [T_ref, T_PI, T_res, error_Wm, Wm_filt] = ...
    PI_speed( ...
        Wm_ref, Wm, Tl_tdm, ...
        params)

% =========================================================================
% CONTROL DE VELOCIDAD PI/PID + CONTROL RESONANTE
%
% Entradas:
%
%   Wm_ref      : referencia de velocidad mecánica [rad/s]
%   Wm          : velocidad mecánica medida/estimada [rad/s]
%   Tl_tdm      : torque de carga/feedforward estimado [Nm]
%
%
% Salidas:
%
%   T_ref       : referencia total de torque [Nm]
%   T_PI        : contribución PI/PID [Nm]
%   T_res       : contribución resonante [Nm]
%   error_Wm    : error de velocidad [rad/s]
%   Wm_filt     : velocidad filtrada [rad/s]
%
%
% Arquitectura:
%
%                     Wm_ref
%                        |
%                        v
%                  error velocidad
%                    /       \
%                   /         \
%              PI/PID       Resonante
%                |              |
%              T_PI           T_res
%                 \             /
%                  \           /
%                   +---------+
%                        |
%              + Tl_tdm + B*Wm
%                        |
%                        v
%                      T_ref
%                        |
%                        v
%                      MTPA
%                        |
%                  Id_ref, Iq_ref
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


% -------------------------------------------------------------------------
% Parámetros resonantes
% -------------------------------------------------------------------------

Kr = params.RI.Kr;

zeta_res = params.RI.zeta;

wr = params.RI.wr;

T_res_max = params.RI.T_res_max;

res_enable = params.RI.enable;

% =========================================================================
% 1. Estados persistentes
% =========================================================================

persistent int_w

persistent Wm_filt_prev
persistent Wm_prev

persistent dWm_filt_prev


% Estados del resonador
persistent xr1
persistent xr2


if isempty(int_w)

    int_w = 0.0;

    Wm_filt_prev = 0.0;
    Wm_prev      = 0.0;

    dWm_filt_prev = 0.0;

    xr1 = 0.0;
    xr2 = 0.0;

end


% =========================================================================
% 2. Protección numérica
% =========================================================================

Ts = max(Ts, 1e-9);

Tf_w = max(Tf_w, Ts);

Tmax = max(abs(Tmax), 1e-9);

wr = max(abs(wr), 1e-6);

zeta_res = max(zeta_res, 1e-6);

T_res_max = max(abs(T_res_max), 1e-9);


% =========================================================================
% 3. Filtro pasabajos de velocidad
% =========================================================================

alpha_w = ...
    Ts / (Tf_w + Ts);


Wm_filt = ...
    Wm_filt_prev ...
    + alpha_w * ...
    (Wm - Wm_filt_prev);


% =========================================================================
% 4. Error de velocidad
% =========================================================================

error_Wm = ...
    Wm_ref - Wm_filt;


% =========================================================================
% 5. Derivada de velocidad
% =========================================================================

dWm_raw = ...
    (Wm_filt - Wm_prev) / Ts;


% =========================================================================
% 6. Filtrado de derivada
% =========================================================================

alpha_d = ...
    Ts / (Tf_w + Ts);


dWm_filt = ...
    dWm_filt_prev ...
    + alpha_d * ...
    (dWm_raw - dWm_filt_prev);


% =========================================================================
% 7. Término proporcional
% =========================================================================

T_P = ...
    Kp_w * error_Wm;


% =========================================================================
% 8. Término derivativo
% =========================================================================

T_D = ...
    -Kd_w * dWm_filt;


% =========================================================================
% 9. Torque producido por PI/PID
% =========================================================================

T_PI = ...
    T_P ...
    + int_w ...
    + T_D;


% =========================================================================
% 10. CONTROL RESONANTE
%
%                    2*zeta*wr*s
% Gres(s) = Kr -------------------------------
%                 s² + 2*zeta*wr*s + wr²
%
%
% Estados:
%
%   xr1_dot = xr2
%
%   xr2_dot =
%
%       -wr²*xr1
%       -2*zeta*wr*xr2
%       +2*zeta*wr*Kr*error_Wm
%
%
%   T_res = xr2
%
% =========================================================================

if res_enable ~= 0

    dxr1 = xr2;

    dxr2 = ...
        -(wr * wr) * xr1 ...
        - 2.0 * zeta_res * wr * xr2 ...
        + 2.0 * zeta_res * wr * Kr * error_Wm;


    % ---------------------------------------------------------------------
    % Integración discreta
    % ---------------------------------------------------------------------

    xr1 = ...
        xr1 ...
        + Ts * dxr1;

    xr2 = ...
        xr2 ...
        + Ts * dxr2;


    % ---------------------------------------------------------------------
    % Salida resonante
    % ---------------------------------------------------------------------

    T_res = xr2;


    % ---------------------------------------------------------------------
    % Saturación independiente del resonador
    % ---------------------------------------------------------------------

    T_res = ...
        min( ...
            max(T_res, -T_res_max), ...
            T_res_max);

else

    % Si se desactiva el resonador,
    % reiniciamos sus estados.

    xr1 = 0.0;
    xr2 = 0.0;

    T_res = 0.0;

end


% =========================================================================
% 11. Feedforward
% =========================================================================

T_ff = ...
    Tl_tdm ...
    + B * Wm_filt;


% =========================================================================
% 12. Referencia de torque antes de saturación
%
% Aquí aparece el cambio importante:
%
%   T_unsat = T_PI + T_res + T_ff
%
% =========================================================================

T_unsat = ...
    T_PI ...
    + T_res ...
    + T_ff;


% =========================================================================
% 13. Saturación de referencia de torque
% =========================================================================

T_ref = ...
    min( ...
        max(T_unsat, -Tmax), ...
        Tmax);


% =========================================================================
% 14. Anti-windup condicional
%
% El anti-windup considera ahora la suma COMPLETA:
%
%   T_PI + T_res + T_ff
%
% porque cualquiera de estos términos puede llevar T_ref a saturación.
%
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
        + Ki_w * error_Wm * Ts;

end


% =========================================================================
% 15. Actualizar estados
% =========================================================================

Wm_filt_prev = Wm_filt;

Wm_prev = Wm_filt;

dWm_filt_prev = dWm_filt;


end