function [Iq_ref, Id_ref, Iq_PI, error_Wm,Wm_filt_out] = ...
    PI_Speed_MTPA(Wm_ref, Wm, params)

% =========================================================================
% CONTROL DE VELOCIDAD Wm
%
% Incluye:
% - PI de velocidad.
% - Derivada opcional filtrada sobre la velocidad medida.
% - Compensación feedforward de carga/fricción.
% - Resonador opcional para cogging.
% - MTPA continuo.
% - Saturación vectorial de corriente.
% - Anti-windup del integrador de velocidad.
%
% Para usar solo PI:
%     params.Kd_w = 0;
%     params.RI_enable = 0;
% =========================================================================


% =========================================================================
% 0. Lectura de parámetros
% =========================================================================
Ts = params.Ts;

Kp_w = params.Kp_w;
Ki_w = params.Ki_w;
Kd_w = params.Kd_w;
Tf_w = params.Tf_w;

Tl_tdm = params.Tl_tdm;
B      = params.B;

Imax = params.Imax;


Ld = params.Ld;
Lq = params.Lq;
Ke = params.Ke;
P  = params.P;


% =========================================================================
% Estados persistentes
% =========================================================================

persistent int_w
persistent Wm_prev dWm_filt_prev
persistent Wm_filt

if isempty(int_w)

    int_w = 0;

    Wm_prev       = 0;
    dWm_filt_prev = 0;
    Wm_filt = 0;

end


% =========================================================================
% Protección numérica
% =========================================================================

Ts_eff   = max(Ts, 1e-9);
Imax_eff = max(abs(Imax), 1e-6);

% =========================================================================
% FILTRO PASABAJOS DE VELOCIDAD
% =========================================================================

fc_Wm = 50;   % Hz, frecuencia de corte del filtro

tau_Wm = 1/(2*pi*fc_Wm);

alpha_Wm = Ts/(Ts + tau_Wm);

Wm_filt = ...
    Wm_filt ...
    + alpha_Wm*(Wm - Wm_filt);

% Constante de torque de reluctancia
Krel = P * (Lq - Ld);


% =========================================================================
% 1. Error de velocidad
% =========================================================================

error_Wm = Wm_ref - Wm_filt;


% =========================================================================
% 2. Derivada filtrada sobre la velocidad medida
%
% Se usa derivada de Wm, no de error, para evitar golpes derivativos
% cuando cambia Wm_ref.
% =========================================================================

dWm_raw = (Wm_filt - Wm_prev) / Ts_eff;

if Tf_w > 0

    alpha_d = Ts_eff / (Tf_w + Ts_eff);

    dWm_filt = ...
        dWm_filt_prev + ...
        alpha_d * (dWm_raw - dWm_filt_prev);

else

    dWm_filt = dWm_raw;

end


% Acción derivativa amortiguadora
D_w = -Kd_w * dWm_filt;


% =========================================================================
% 3. Feedforward de torque de carga y fricción
% =========================================================================

T_ff = Tl_tdm + B * Wm_ref;

if abs(Ke) > 1e-9

    Iq_ff = T_ff / Ke;

else

    Iq_ff = 0;

end


% =========================================================================
% 4. Resonador para compensación de ripple/cogging
% =========================================================================



% =========================================================================
% 5. Integrador candidato
% =========================================================================

int_w_candidate = ...
    int_w + Ki_w * Ts_eff * error_Wm;


% =========================================================================
% 6. Cálculo PI/PID, MTPA y límite de corriente
%
% Se hacen hasta dos pasadas:
%
% Primera:
%   prueba integrando.
%
% Segunda:
%   si saturó y el integrador empeora la saturación,
%   se rechaza esa integración.
% =========================================================================

for pass = 1:2

    % ---------------------------------------------------------------------
    % Salida del controlador de velocidad
    % ---------------------------------------------------------------------

    Iq_PI = ...
        Kp_w * error_Wm ...
        + int_w_candidate ...
        + D_w;


    % ---------------------------------------------------------------------
    % Corriente q total solicitada
    % ---------------------------------------------------------------------

    Iq_ref_raw = ...
        Iq_ff ...
        + Iq_PI;


    % =====================================================================
    % MTPA continuo
    %
    % Id_ref tiende suavemente a cero cuando Iq_ref_raw tiende a cero.
    % No existe umbral que fuerce saltos Id_ref = 0.
    % =====================================================================

    if abs(Krel) > 1e-9 && abs(Ke) > 1e-9

        A_mtpa = Ke / (2 * Krel);

   

            % Forma numéricamente estable de:
            %
            % Id = A - sign(A)*sqrt(A^2 + Iq^2)
            %

            Id_ref_raw = ...
    -sign(A_mtpa) * (Iq_ref_raw^2) / ...
    (abs(A_mtpa) + sqrt(A_mtpa^2 + Iq_ref_raw^2));

        

    else

        Id_ref_raw = 0;

    end


    % =====================================================================
    % Saturación vectorial de corriente
    % =====================================================================

    Is_mag = sqrt( ...
        Id_ref_raw^2 + ...
        Iq_ref_raw^2);

    if Is_mag > Imax_eff

        ratio = Imax_eff / Is_mag;

        Id_ref = Id_ref_raw * ratio;
        Iq_ref = Iq_ref_raw * ratio;

        sat = true;

    else

        Id_ref = Id_ref_raw;
        Iq_ref = Iq_ref_raw;

        sat = false;

    end


    % =====================================================================
    % Anti-windup condicional
    %
    % Si la salida está saturada y el error pide aumentar todavía más
    % la corriente saturada, no se integra.
    %
    % Si el error ayuda a salir de saturación,
    % sí se permite integrar.
    % =====================================================================

    integrator_worsens_sat = ...
        (Iq_ref_raw * error_Wm) >= 0;

    if pass == 1 && sat && integrator_worsens_sat

        % Rechaza la nueva integración
        % y recalcula la salida.

        int_w_candidate = int_w;

    else

        break;

    end

end


% =========================================================================
% Guarda el integrador final
% =========================================================================

int_w = int_w_candidate;

% =========================================================================
% 7. Actualización de memorias
% =========================================================================

Wm_prev         = Wm_filt;
dWm_filt_prev   = dWm_filt;
Wm_filt_out = Wm_filt;


end