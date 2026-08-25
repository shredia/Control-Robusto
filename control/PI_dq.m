function [Vd, Vq, Id_out, Iq_out, ...
          Valpha, Vbetha, Vd_hfi, Vq_hfi, ...
          error_d, error_q, ...
          Id_ref_out, Iq_ref_out, T_ref_real] = ...
    PI_dq( ...
        Wm, X, Theta_e, ...
        T_ref, ...
        Iq_comp,...
        params)

% =========================================================================
% PI dq VECTORIAL + MTPA + INYECCIÓN HFI EN EJE d
%
% Entradas variables:
%
%   Wm          : velocidad mecánica [rad/s]
%   X           : [Ia; Ib; ...]
%   Theta_e     : ángulo eléctrico utilizado por el controlador [rad]
%   T_ref       : referencia de torque entregada por control de velocidad [Nm]
%
%
% Parámetros:
%
%   params.Ts
%   params.Vdc
%
%   params.Kp_d
%   params.Ki_d
%   params.Kp_q
%   params.Ki_q
%
%   params.Psi
%   params.Ld
%   params.Lq
%   params.P
%
%   params.Imax
%
%   params.Amplitud_HFI
%   params.wh
%   params.HFI_flag
%
%
% Salidas:
%
%   Vd, Vq
%       Tensión fundamental dq saturada, SIN HFI
%
%   Vd_hfi, Vq_hfi
%       Tensión dq total incluyendo HFI
%
%   Valpha, Vbetha
%       Tensión alpha-beta total
%
%   Id_out, Iq_out
%       Corrientes medidas en marco dq
%
%   error_d, error_q
%       Errores de los PI
%
%   Id_ref_out, Iq_ref_out
%       Referencias calculadas internamente mediante MTPA
%
%   T_ref_real
%       Torque correspondiente a Id_ref/Iq_ref después de saturar Imax
%
%
% =========================================================================
%
% MODELO DE TORQUE UTILIZADO
%
%   Te = P*(Psi + (Ld-Lq)*Id)*Iq
%
% Como:
%
%   DeltaL = Lq - Ld
%
% entonces:
%
%   Te = P*(Psi - DeltaL*Id)*Iq
%
%
% Para el motor actual:
%
%   Lq > Ld
%
% por lo tanto:
%
%   DeltaL > 0
%
% y para torque positivo el MTPA normalmente produce:
%
%   Id_ref < 0
%   Iq_ref > 0
%
% =========================================================================


% =========================================================================
% 0. Parámetros
% =========================================================================

Ts = params.Ts;

Vdc = params.Vdc;

Kp_d = params.Kp_d_salient;
Ki_d = params.Ki_d_salient;

Kp_q = params.Kp_q_salient;
Ki_q = params.Ki_q_salient;

Psi = params.Psi;

Ld = params.Ld;
Lq = params.Lq;

P = params.P;

Imax = params.Imax;

Amplitud_HFI = params.Amplitud_HFI;
wh           = params.wh;
HFI_enable   = params.HFI_flag;


% =========================================================================
% Estados persistentes
% =========================================================================

persistent v_int_d v_int_q
persistent t_hfi

if isempty(v_int_d)

    v_int_d = 0.0;
    v_int_q = 0.0;

    t_hfi = 0.0;

end


% =========================================================================
% Protección básica
% =========================================================================

Imax = max(abs(Imax), 1e-6);
Ts   = max(Ts, 1e-9);

P_eff = max(abs(P), 1e-9);


% =========================================================================
% 1. Transformada de Park
% =========================================================================

Ia = X(1);
Ib = X(2);

cos_th = cos(Theta_e);
sin_th = sin(Theta_e);

Id =  Ia*cos_th + Ib*sin_th;
Iq = -Ia*sin_th + Ib*cos_th;


% =========================================================================
% 2. MTPA
%
% Se busca:
%
%       minimizar:
%
%           Is^2 = Id^2 + Iq^2
%
%       sujeto a:
%
%           T_ref = P*(Psi - DeltaL*Id)*Iq
%
%
% donde:
%
%       DeltaL = Lq - Ld
%
%
% La condición MTPA resulta:
%
%       DeltaL*Iq^2
%       + Psi*Id
%       - DeltaL*Id^2 = 0
%
%
% Se resuelve numéricamente usando Newton-Raphson sobre Id.
%
% =========================================================================

DeltaL = Lq - Ld;


% -------------------------------------------------------------------------
% Torque prácticamente nulo
% -------------------------------------------------------------------------

if abs(T_ref) < 1e-10

    Id_ref = 0.0;
    Iq_ref = 0.0;


% -------------------------------------------------------------------------
% Motor prácticamente no saliente
%
% Ld ~= Lq
%
% En ese caso MTPA:
%
% Id_ref = 0
% Iq_ref = T/(P*Psi)
% -------------------------------------------------------------------------

elseif abs(DeltaL) < 1e-10

    Id_ref = 0.0;

    torque_constant = P_eff * Psi;

    if abs(torque_constant) < 1e-9

        Iq_ref = 0.0;

    else

        Iq_ref = T_ref / torque_constant;

    end


% -------------------------------------------------------------------------
% Motor saliente
% -------------------------------------------------------------------------

else

    % =====================================================================
    % 2.1 Estimación inicial
    %
    % Ignorando inicialmente el torque de reluctancia:
    %
    %       Iq0 = T/(P*Psi)
    %
    % =====================================================================

    torque_constant = P_eff * Psi;

    if abs(torque_constant) > 1e-9

        Iq0 = T_ref / torque_constant;

    else

        Iq0 = 0.0;

    end


    % =====================================================================
    % Estimación analítica inicial de Id
    %
    % De la condición:
    %
    % DeltaL*Iq^2 + Psi*Id - DeltaL*Id^2 = 0
    %
    % para Iq = Iq0.
    %
    % La raíz físicamente interesante es la que tiende a Id=0
    % cuando Iq -> 0.
    % =====================================================================

    ratio = Psi / DeltaL;

    Id_ref = ...
        0.5 * ...
        (ratio - sign(ratio) * ...
        sqrt(ratio*ratio + 4.0*Iq0*Iq0));


    % Protección
    Id_ref = min(max(Id_ref, -Imax), Imax);


    % =====================================================================
    % 2.2 Newton-Raphson
    %
    % Torque:
    %
    % Iq(Id) =
    %
    %       T_ref
    % -------------------
    % P*(Psi-DeltaL*Id)
    %
    %
    % Condición MTPA:
    %
    % f(Id) =
    %
    % DeltaL*Iq^2
    % + Psi*Id
    % - DeltaL*Id^2
    %
    % =====================================================================

    for k = 1:6

        flux_eff = ...
            Psi - DeltaL*Id_ref;


        % -------------------------------------------------------------
        % Evitar división por cero
        % -------------------------------------------------------------

        if abs(flux_eff) < 1e-8

            if flux_eff >= 0.0

                flux_eff = 1e-8;

            else

                flux_eff = -1e-8;

            end

        end


        % -------------------------------------------------------------
        % Iq correspondiente al torque solicitado
        % -------------------------------------------------------------

        Iq_tmp = ...
            T_ref / ...
            (P_eff * flux_eff);


        % -------------------------------------------------------------
        % Función MTPA
        % -------------------------------------------------------------

        f = ...
            DeltaL * Iq_tmp * Iq_tmp ...
            + Psi * Id_ref ...
            - DeltaL * Id_ref * Id_ref;


        % -------------------------------------------------------------
        % dIq / dId
        %
        % Iq =
        %
        % T/[P*(Psi-DeltaL*Id)]
        %
        % -------------------------------------------------------------

        dIq_dId = ...
            T_ref * DeltaL / ...
            (P_eff * flux_eff * flux_eff);


        % -------------------------------------------------------------
        % df / dId
        % -------------------------------------------------------------

        df = ...
            2.0 * DeltaL * Iq_tmp * dIq_dId ...
            + Psi ...
            - 2.0 * DeltaL * Id_ref;


        % -------------------------------------------------------------
        % Newton
        % -------------------------------------------------------------

        if abs(df) > 1e-10

            Id_new = Id_ref - f/df;

        else

            Id_new = Id_ref;

        end


        % -------------------------------------------------------------
        % Protección numérica
        % -------------------------------------------------------------

        Id_new = ...
            min(max(Id_new, -Imax), Imax);


        Id_ref = Id_new;

    end


    % =====================================================================
    % 2.3 Calcular Iq_ref final usando la ecuación de torque
    % =====================================================================

    flux_eff = ...
        Psi - DeltaL*Id_ref;


    if abs(flux_eff) < 1e-8

        if flux_eff >= 0.0

            flux_eff = 1e-8;

        else

            flux_eff = -1e-8;

        end

    end


    Iq_ref = ...
        T_ref / ...
        (P_eff * flux_eff);

end


% =========================================================================
% 3. Saturación vectorial de referencias de corriente
%
% Se conserva la dirección del vector MTPA:
%
%        sqrt(Id_ref^2 + Iq_ref^2) <= Imax
%
% =========================================================================


Iq_ref = ...
    Iq_ref + Iq_comp;

I_ref_mag_raw = hypot(Id_ref, Iq_ref);


if I_ref_mag_raw > Imax

    scale_I = ...
        Imax/I_ref_mag_raw;

    Id_ref_lim = ...
        Id_ref * scale_I;

    Iq_ref_lim = ...
        Iq_ref * scale_I;

else

    Id_ref_lim = Id_ref;
    Iq_ref_lim = Iq_ref;

end


% =========================================================================
% 4. Torque realmente alcanzable después del límite de corriente
% =========================================================================

T_ref_real = ...
    P * ...
    (Psi - DeltaL*Id_ref_lim) * ...
    Iq_ref_lim;


% =========================================================================
% 5. Salidas de referencia MTPA
% =========================================================================

Id_ref_out = Id_ref_lim;
Iq_ref_out = Iq_ref_lim;


% =========================================================================
% 6. Errores de corriente
% =========================================================================

err_d = ...
    Id_ref_lim - Id;

err_q = ...
    Iq_ref_lim - Iq;


error_d = err_d;
error_q = err_q;


% =========================================================================
% 7. Velocidad eléctrica
% =========================================================================

We = ...
    P * Wm;


% =========================================================================
% 8. Desacople y feedforward
%
% Modelo:
%
% Vd =
%
% R*Id
% + Ld*dId/dt
% - We*Lq*Iq
%
%
% Vq =
%
% R*Iq
% + Lq*dIq/dt
% + We*Ld*Id
% + We*Psi
%
%
% Se utilizan las referencias limitadas para el desacople.
% =========================================================================

V_cross_d = ...
    -We * Lq * Iq_ref_lim;


V_cross_q = ...
     We * Ld * Id_ref_lim ...
   + We * Psi;


% =========================================================================
% 9. PI sin saturar
% =========================================================================

Vd_unsat = ...
    Kp_d * err_d ...
    + v_int_d ...
    + V_cross_d;


Vq_unsat = ...
    Kp_q * err_q ...
    + v_int_q ...
    + V_cross_q;


% =========================================================================
% 10. Límite de tensión disponible
% =========================================================================

% Actualmente:
%
%       Vmax = Vdc
%
% Dependiendo de tu etapa de potencia/modulación puede posteriormente
% corresponder:
%
%       Vdc/sqrt(2)
%
%       Vdc/sqrt(3)
%
% etc.
%
% Mantengo el comportamiento de tu controlador actual.

Vmax_inv = Vdc;


% =========================================================================
% Reserva de tensión para HFI
% =========================================================================

if HFI_enable == 1

    Vmax_control = ...
        max( ...
        Vmax_inv - abs(Amplitud_HFI), ...
        0.0);

else

    Vmax_control = Vmax_inv;

end


% =========================================================================
% 11. Saturación vectorial de tensión fundamental
% =========================================================================

V_mag_unsat = ...
    hypot(Vd_unsat, Vq_unsat);


if V_mag_unsat > Vmax_control && ...
        V_mag_unsat > 1e-9

    scale_V = ...
        Vmax_control / V_mag_unsat;


    Vd = ...
        Vd_unsat * scale_V;

    Vq = ...
        Vq_unsat * scale_V;

else

    Vd = Vd_unsat;
    Vq = Vq_unsat;

end


% =========================================================================
% 12. Anti-windup por back-calculation
% =========================================================================

Kaw_d = ...
    Ki_d / max(abs(Kp_d), 1e-9);

Kaw_q = ...
    Ki_q / max(abs(Kp_q), 1e-9);


v_int_d = ...
    v_int_d ...
    + Ts * ( ...
        Ki_d * err_d ...
        + Kaw_d * (Vd - Vd_unsat));


v_int_q = ...
    v_int_q ...
    + Ts * ( ...
        Ki_q * err_q ...
        + Kaw_q * (Vq - Vq_unsat));


% =========================================================================
% 13. Inyección HFI
% =========================================================================

t_hfi = ...
    t_hfi + Ts;


% =========================================================================
% Evitar crecimiento indefinido de t_hfi
% =========================================================================

if wh > 0.0

    T_hfi = ...
        2.0*pi / wh;


    if t_hfi >= T_hfi

        t_hfi = ...
            t_hfi - T_hfi;

    end

end


% =========================================================================
% Señal HFI
% =========================================================================

if HFI_enable == 1

    v_hfi = ...
        Amplitud_HFI * ...
        sin(wh*t_hfi);

else

    v_hfi = 0.0;

end


% =========================================================================
% Inyección HFI sobre eje d estimado
% =========================================================================

Vd_hfi = ...
    Vd + v_hfi;

Vq_hfi = ...
    Vq;


% =========================================================================
% 14. Saturación final incluyendo HFI
% =========================================================================

V_mag_hfi = ...
    hypot(Vd_hfi, Vq_hfi);


if V_mag_hfi > Vmax_inv && ...
        V_mag_hfi > 1e-9

    scale_hfi = ...
        Vmax_inv / V_mag_hfi;


    Vd_hfi = ...
        Vd_hfi * scale_hfi;

    Vq_hfi = ...
        Vq_hfi * scale_hfi;

end


% =========================================================================
% 15. Transformada inversa de Park
% =========================================================================

Valpha = ...
    Vd_hfi*cos_th ...
    - Vq_hfi*sin_th;


Vbetha = ...
    Vd_hfi*sin_th ...
    + Vq_hfi*cos_th;


% =========================================================================
% 16. Salidas de monitoreo
% =========================================================================

Id_out = Id;
Iq_out = Iq;

end