function [Vd, Vq, Id_out, Iq_out, ...
          Valpha, Vbetha, Vd_hfi, Vq_hfi, ...
          error_d, error_q, ...
          Id_ref_out, Iq_ref_out, T_ref_real] = ...
    PI_dq( ...
        Wm, X, Theta_e, ...
        T_ref, ...
        comp_res, ...
        params)

% =========================================================================
% PI dq VECTORIAL + MTPA + PR d/q + HFI
%
% -------------------------------------------------------------------------
% CONTROL DE CORRIENTE
%
% Eje d:
%
%   Vd =
%       PI_d(error_d)
%       + PR_d(error_d)
%       + desacople_d
%
%
% Eje q:
%
%   Vq =
%       PI_q(error_q)
%       + PR_q(error_q)
%       + desacople_q
%       + comp_res
%
%
% Resonador:
%
%                       2*Kr*wc*s
%       Gres(s) = -----------------------
%                  s^2 + 2*wc*s + wr^2
%
%
% Frecuencia resonante:
%
%       wr = Nc * abs(Wm)
%
% donde Nc corresponde al número de períodos de cogging por revolución
% mecánica.
%
%
% Discretización:
%
%       Tustin / transformación bilineal
%
%                  2       1-z^-1
%       s = ------------- ----------
%                  Ts      1+z^-1
%
%
% =========================================================================


% =========================================================================
% 0. PARÁMETROS
% =========================================================================

Ts = params.Ts;

Vdc = params.Vdc;


% -------------------------------------------------------------------------
% PI corriente
% -------------------------------------------------------------------------

Kp_d = params.Kp_d_salient;
Ki_d = params.Ki_d_salient;

Kp_q = params.Kp_q_salient;
Ki_q = params.Ki_q_salient;


% -------------------------------------------------------------------------
% Motor
% -------------------------------------------------------------------------

Psi = params.Psi;

Ld = params.Ld;
Lq = params.Lq;

P = params.P;

Imax = params.Imax;


% -------------------------------------------------------------------------
% HFI
% -------------------------------------------------------------------------

Amplitud_HFI = params.Amplitud_HFI;

wh = params.wh;

HFI_enable = params.HFI_flag;


% -------------------------------------------------------------------------
% PR de corriente
% -------------------------------------------------------------------------

PR_enable = params.PR_i.enable;

Nc = params.PR_i.Nc;

wc = params.PR_i.wc;

Kr_d = params.PR_i.Kr_d;
Kr_q = params.PR_i.Kr_q;

Vres_max_d = params.PR_i.Vres_max_d;
Vres_max_q = params.PR_i.Vres_max_q;

wr_min = params.PR_i.wr_min;
wr_max = params.PR_i.wr_max;


% =========================================================================
% 1. ESTADOS PERSISTENTES
% =========================================================================

persistent v_int_d
persistent v_int_q

persistent t_hfi


% -------------------------------------------------------------------------
% Estados PR eje d
%
% e_d[k-1]
% e_d[k-2]
%
% y_res_d[k-1]
% y_res_d[k-2]
% -------------------------------------------------------------------------

persistent ed_z1
persistent ed_z2

persistent vdres_z1
persistent vdres_z2


% -------------------------------------------------------------------------
% Estados PR eje q
% -------------------------------------------------------------------------

persistent eq_z1
persistent eq_z2

persistent vqres_z1
persistent vqres_z2


% =========================================================================
% Inicialización
% =========================================================================

if isempty(v_int_d)

    % PI
    v_int_d = 0.0;
    v_int_q = 0.0;

    % HFI
    t_hfi = 0.0;

    % PR d
    ed_z1 = 0.0;
    ed_z2 = 0.0;

    vdres_z1 = 0.0;
    vdres_z2 = 0.0;

    % PR q
    eq_z1 = 0.0;
    eq_z2 = 0.0;

    vqres_z1 = 0.0;
    vqres_z2 = 0.0;

end


% =========================================================================
% 2. PROTECCIÓN NUMÉRICA
% =========================================================================

Imax = ...
    max(abs(Imax), 1e-6);

Ts = ...
    max(Ts, 1e-9);

P_eff = ...
    max(abs(P), 1e-9);

wc = ...
    max(abs(wc), 1e-6);

wr_min = ...
    max(abs(wr_min), 0.0);

wr_max = ...
    max(abs(wr_max), wr_min);


% =========================================================================
% 3. TRANSFORMADA DE PARK
% =========================================================================

Ia = X(1);
Ib = X(2);


cos_th = cos(Theta_e);
sin_th = sin(Theta_e);


Id = ...
     Ia*cos_th ...
   + Ib*sin_th;


Iq = ...
    -Ia*sin_th ...
    + Ib*cos_th;


% =========================================================================
% 4. MTPA
%
% Modelo de torque:
%
%       Te = P*(Psi + (Ld-Lq)*Id)*Iq
%
% Como:
%
%       DeltaL = Lq - Ld
%
% entonces:
%
%       Te = P*(Psi - DeltaL*Id)*Iq
%
%
% MTPA:
%
%       minimizar:
%
%           Id^2 + Iq^2
%
%       sujeto a:
%
%           Tref = P*(Psi - DeltaL*Id)*Iq
%
%
% Condición:
%
%       DeltaL*Iq^2
%       + Psi*Id
%       - DeltaL*Id^2 = 0
%
% =========================================================================

DeltaL = ...
    Lq - Ld;


% -------------------------------------------------------------------------
% Torque prácticamente nulo
% -------------------------------------------------------------------------

if abs(T_ref) < 1e-10

    Id_ref = 0.0;
    Iq_ref = 0.0;


% -------------------------------------------------------------------------
% Motor prácticamente no saliente
% -------------------------------------------------------------------------

elseif abs(DeltaL) < 1e-10

    Id_ref = 0.0;

    torque_constant = ...
        P_eff * Psi;


    if abs(torque_constant) < 1e-9

        Iq_ref = 0.0;

    else

        Iq_ref = ...
            T_ref / torque_constant;

    end


% -------------------------------------------------------------------------
% Motor saliente
% -------------------------------------------------------------------------

else

    % =====================================================================
    % 4.1 Estimación inicial de Iq
    % =====================================================================

    torque_constant = ...
        P_eff * Psi;


    if abs(torque_constant) > 1e-9

        Iq0 = ...
            T_ref / torque_constant;

    else

        Iq0 = 0.0;

    end


    % =====================================================================
    % 4.2 Estimación analítica inicial de Id
    % =====================================================================

    ratio = ...
        Psi / DeltaL;


    Id_ref = ...
        0.5 * ...
        ( ...
        ratio ...
        - sign(ratio) * ...
        sqrt( ...
        ratio*ratio ...
        + 4.0*Iq0*Iq0) ...
        );


    % Protección
    Id_ref = ...
        min( ...
        max(Id_ref, -Imax), ...
        Imax);


    % =====================================================================
    % 4.3 Newton-Raphson
    % =====================================================================

    for k = 1:6

        flux_eff = ...
            Psi ...
            - DeltaL*Id_ref;


        % -----------------------------------------------------------------
        % Evitar división por cero
        % -----------------------------------------------------------------

        if abs(flux_eff) < 1e-8

            if flux_eff >= 0.0

                flux_eff = 1e-8;

            else

                flux_eff = -1e-8;

            end

        end


        % -----------------------------------------------------------------
        % Iq correspondiente al torque solicitado
        % -----------------------------------------------------------------

        Iq_tmp = ...
            T_ref / ...
            (P_eff * flux_eff);


        % -----------------------------------------------------------------
        % Función MTPA
        % -----------------------------------------------------------------

        f = ...
            DeltaL * Iq_tmp * Iq_tmp ...
            + Psi * Id_ref ...
            - DeltaL * Id_ref * Id_ref;


        % -----------------------------------------------------------------
        % dIq / dId
        % -----------------------------------------------------------------

        dIq_dId = ...
            T_ref * DeltaL / ...
            ( ...
            P_eff ...
            * flux_eff ...
            * flux_eff);


        % -----------------------------------------------------------------
        % df / dId
        % -----------------------------------------------------------------

        df = ...
            2.0 * DeltaL ...
            * Iq_tmp ...
            * dIq_dId ...
            + Psi ...
            - 2.0 * DeltaL * Id_ref;


        % -----------------------------------------------------------------
        % Newton
        % -----------------------------------------------------------------

        if abs(df) > 1e-10

            Id_new = ...
                Id_ref - f/df;

        else

            Id_new = ...
                Id_ref;

        end


        % -----------------------------------------------------------------
        % Protección
        % -----------------------------------------------------------------

        Id_new = ...
            min( ...
            max(Id_new, -Imax), ...
            Imax);


        Id_ref = ...
            Id_new;

    end


    % =====================================================================
    % 4.4 Iq final a partir del torque
    % =====================================================================

    flux_eff = ...
        Psi ...
        - DeltaL*Id_ref;


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
% 5. SATURACIÓN VECTORIAL DE REFERENCIAS DE CORRIENTE
% =========================================================================

I_ref_mag_raw = ...
    hypot(Id_ref, Iq_ref);


if I_ref_mag_raw > Imax

    scale_I = ...
        Imax / I_ref_mag_raw;


    Id_ref_lim = ...
        Id_ref * scale_I;


    Iq_ref_lim = ...
        Iq_ref * scale_I;

else

    Id_ref_lim = ...
        Id_ref;


    Iq_ref_lim = ...
        Iq_ref;

end


% =========================================================================
% 6. TORQUE REALMENTE SOLICITABLE
% =========================================================================

T_ref_real = ...
    P ...
    * (Psi - DeltaL*Id_ref_lim) ...
    * Iq_ref_lim;


% =========================================================================
% 7. SALIDAS REFERENCIA MTPA
% =========================================================================

Id_ref_out = ...
    Id_ref_lim;


Iq_ref_out = ...
    Iq_ref_lim;


% =========================================================================
% 8. ERRORES DE CORRIENTE
% =========================================================================

err_d = ...
    Id_ref_lim - Id;


err_q = ...
    Iq_ref_lim - Iq;


error_d = ...
    err_d;


error_q = ...
    err_q;


% =========================================================================
% 9. FRECUENCIA RESONANTE
%
% Cogging:
%
%       wr = Nc * |Wm|
%
% Ejemplo:
%
%       Nc = 200
%       Wm = 10 rad/s
%
%       wr = 2000 rad/s
%       fr = 318.3 Hz
%
% =========================================================================

wr = ...
    Nc * abs(Wm);


% -------------------------------------------------------------------------
% Limitar frecuencia utilizada por los PR
% -------------------------------------------------------------------------

wr = ...
    min( ...
    max(wr, wr_min), ...
    wr_max);


% =========================================================================
% 10. CONTROL RESONANTE DE CORRIENTE
%
%                         2*Kr*wc*s
%       Gres(s) = -------------------------
%                   s² + 2*wc*s + wr²
%
%
% Tustin:
%
%              2    1-z^-1
%       s = ------- --------
%              Ts   1+z^-1
%
%
% Después de sustituir:
%
%       y[k] =
%
%       b0*e[k]
%       + b2*e[k-2]
%       - a1*y[k-1]
%       - a2*y[k-2]
%       --------------------
%              a0
%
% donde:
%
%       b0 =  2*Kr*wc*K
%       b1 =  0
%       b2 = -2*Kr*wc*K
%
%       K = 2/Ts
%
% =========================================================================

Vd_res = 0.0;
Vq_res = 0.0;


if PR_enable == 1

    % =====================================================================
    % Constante Tustin
    % =====================================================================

    Kt = ...
        2.0 / Ts;


    % =====================================================================
    % Denominador
    % =====================================================================

    a0 = ...
        Kt*Kt ...
        + 2.0*wc*Kt ...
        + wr*wr;


    a1 = ...
        -2.0*Kt*Kt ...
        + 2.0*wr*wr;


    a2 = ...
        Kt*Kt ...
        - 2.0*wc*Kt ...
        + wr*wr;


    % =====================================================================
    % PR EJE d
    % =====================================================================

    b0_d = ...
        2.0 ...
        * Kr_d ...
        * wc ...
        * Kt;


    b2_d = ...
        -b0_d;


    % ---------------------------------------------------------------------
    % Salida resonante sin saturar
    % ---------------------------------------------------------------------

    Vd_res_raw = ...
        ( ...
          b0_d * err_d ...
        + b2_d * ed_z2 ...
        - a1   * vdres_z1 ...
        - a2   * vdres_z2 ...
        ) ...
        / a0;


    % ---------------------------------------------------------------------
    % Protección individual
    % ---------------------------------------------------------------------

    Vd_res = ...
        min( ...
        max(Vd_res_raw, -Vres_max_d), ...
        Vres_max_d);


    % =====================================================================
    % PR EJE q
    % =====================================================================

    b0_q = ...
        2.0 ...
        * Kr_q ...
        * wc ...
        * Kt;


    b2_q = ...
        -b0_q;


    Vq_res_raw = ...
        ( ...
          b0_q * err_q ...
        + b2_q * eq_z2 ...
        - a1   * vqres_z1 ...
        - a2   * vqres_z2 ...
        ) ...
        / a0;


    Vq_res = ...
        min( ...
        max(Vq_res_raw, -Vres_max_q), ...
        Vres_max_q);


    % =====================================================================
    % ACTUALIZACIÓN ESTADOS PR
    %
    % Se guarda la salida limitada para evitar crecimiento interno
    % indefinido cuando el resonador alcanza su límite.
    % =====================================================================

    ed_z2 = ...
        ed_z1;

    ed_z1 = ...
        err_d;


    vdres_z2 = ...
        vdres_z1;

    vdres_z1 = ...
        Vd_res;


    eq_z2 = ...
        eq_z1;

    eq_z1 = ...
        err_q;


    vqres_z2 = ...
        vqres_z1;

    vqres_z1 = ...
        Vq_res;


else

    % =====================================================================
    % PR deshabilitado
    % =====================================================================

    Vd_res = 0.0;
    Vq_res = 0.0;


    % Reiniciar estados
    ed_z1 = 0.0;
    ed_z2 = 0.0;

    vdres_z1 = 0.0;
    vdres_z2 = 0.0;


    eq_z1 = 0.0;
    eq_z2 = 0.0;

    vqres_z1 = 0.0;
    vqres_z2 = 0.0;

end


% =========================================================================
% 11. VELOCIDAD ELÉCTRICA
% =========================================================================

We = ...
    P * Wm;


% =========================================================================
% 12. DESACOPLE / FEEDFORWARD
%
% Modelo:
%
%       Vd =
%           R*Id
%           + Ld*dId/dt
%           - We*Lq*Iq
%
%
%       Vq =
%           R*Iq
%           + Lq*dIq/dt
%           + We*Ld*Id
%           + We*Psi
%
%
% Usamos las referencias limitadas.
%
% =========================================================================

V_cross_d = ...
    -We ...
    * Lq ...
    * Iq_ref_lim;


V_cross_q = ...
      We ...
    * Ld ...
    * Id_ref_lim ...
    + We ...
    * Psi;


% =========================================================================
% 13. PI + PR + DESACOPLE
%
%
% Eje d:
%
%       Vd = PI_d + PR_d + desacople
%
%
% Eje q:
%
%       Vq = PI_q + PR_q + desacople + comp_res
%
% =========================================================================

Vd_unsat = ...
      Kp_d * err_d ...
    + v_int_d ...
    + Vd_res ...
    + V_cross_d;


Vq_unsat = ...
      Kp_q * err_q ...
    + v_int_q ...
    + Vq_res ...
    + V_cross_q ...
    + comp_res;


% =========================================================================
% 14. LÍMITE DE TENSIÓN DEL INVERSOR
%
% Actualmente se conserva:
%
%       Vmax = Vdc
%
% Si tu modulación requiere:
%
%       Vdc/sqrt(2)
%
% o:
%
%       Vdc/sqrt(3)
%
% habría que modificar esta expresión.
%
% =========================================================================

Vmax_inv = ...
    Vdc;


% =========================================================================
% 15. RESERVA DE TENSIÓN PARA HFI
% =========================================================================

if HFI_enable == 1

    Vmax_control = ...
        max( ...
        Vmax_inv ...
        - abs(Amplitud_HFI), ...
        0.0);

else

    Vmax_control = ...
        Vmax_inv;

end


% =========================================================================
% 16. SATURACIÓN VECTORIAL DE TENSIÓN FUNDAMENTAL
% =========================================================================

V_mag_unsat = ...
    hypot(Vd_unsat, Vq_unsat);


if V_mag_unsat > Vmax_control ...
        && V_mag_unsat > 1e-9

    scale_V = ...
        Vmax_control ...
        / V_mag_unsat;


    Vd = ...
        Vd_unsat ...
        * scale_V;


    Vq = ...
        Vq_unsat ...
        * scale_V;

else

    Vd = ...
        Vd_unsat;


    Vq = ...
        Vq_unsat;

end


% =========================================================================
% 17. ANTI-WINDUP PI
%
% El anti-windup se mantiene sobre los integradores PI.
%
% El PR tiene su propia limitación:
%
%       +/- Vres_max
%
% =========================================================================

Kaw_d = ...
    Ki_d ...
    / max(abs(Kp_d), 1e-9);


Kaw_q = ...
    Ki_q ...
    / max(abs(Kp_q), 1e-9);


v_int_d = ...
    v_int_d ...
    + Ts * ...
    ( ...
      Ki_d * err_d ...
    + Kaw_d * ...
      (Vd - Vd_unsat) ...
    );


v_int_q = ...
    v_int_q ...
    + Ts * ...
    ( ...
      Ki_q * err_q ...
    + Kaw_q * ...
      (Vq - Vq_unsat) ...
    );


% =========================================================================
% 18. HFI
% =========================================================================

t_hfi = ...
    t_hfi + Ts;


% -------------------------------------------------------------------------
% Evitar crecimiento indefinido
% -------------------------------------------------------------------------

if wh > 0.0

    T_hfi = ...
        2.0*pi / wh;


    if t_hfi >= T_hfi

        t_hfi = ...
            t_hfi - T_hfi;

    end

end


% =========================================================================
% 19. SEÑAL HFI
% =========================================================================

if HFI_enable == 1

    v_hfi = ...
        Amplitud_HFI ...
        * sin(wh*t_hfi);

else

    v_hfi = 0.0;

end


% =========================================================================
% 20. HFI EN EJE d
% =========================================================================

Vd_hfi = ...
    Vd + v_hfi;


Vq_hfi = ...
    Vq;


% =========================================================================
% 21. SATURACIÓN FINAL CON HFI
% =========================================================================

V_mag_hfi = ...
    hypot(Vd_hfi, Vq_hfi);


if V_mag_hfi > Vmax_inv ...
        && V_mag_hfi > 1e-9

    scale_hfi = ...
        Vmax_inv ...
        / V_mag_hfi;


    Vd_hfi = ...
        Vd_hfi ...
        * scale_hfi;


    Vq_hfi = ...
        Vq_hfi ...
        * scale_hfi;

end


% =========================================================================
% 22. TRANSFORMADA INVERSA DE PARK
% =========================================================================

Valpha = ...
      Vd_hfi*cos_th ...
    - Vq_hfi*sin_th;


Vbetha = ...
      Vd_hfi*sin_th ...
    + Vq_hfi*cos_th;


% =========================================================================
% 23. SALIDAS DE MONITOREO
% =========================================================================

Id_out = ...
    Id;


Iq_out = ...
    Iq;


end