function [Vd, Vq, Id_out, Iq_out, ...
          Valpha, Vbetha, Vd_hfi, Vq_hfi, ...
          error_d, error_q] = ...
    PI_dq( ...
        Wm, X, Theta_e, ...
        Iq_ref_req, Id_ref, Iq_comp, ...
        params)

% =========================================================================
% PI dq VECTORIAL + INYECCIÓN HFI EN EJE d
%
% Entradas variables:
%   Wm          : velocidad mecánica
%   X           : [Ia; Ib; ...]
%   Theta_e     : ángulo eléctrico utilizado por el controlador
%   Iq_ref_req  : referencia de corriente q
%   Id_ref      : referencia de corriente d
%   Iq_comp     : compensación adicional de corriente q
%
% Parámetros:
%   params.PI...
%
% Salidas:
%   Vd, Vq             : tensión fundamental saturada, sin HFI
%   Vd_hfi, Vq_hfi     : tensión dq total, con HFI
%   Valpha, Vbetha     : tensión alpha-beta con HFI
%   Id_out, Iq_out     : corrientes medidas en marco dq estimado
%   error_d, error_q   : errores utilizados por los PI
% =========================================================================


% =========================================================================
% 0. Parámetros
% =========================================================================

Ts = params.Ts;

Vdc = params.Vdc;

Kp_d = params.Kp_d;
Ki_d = params.Ki_d;

Kp_q = params.Kp_q;
Ki_q = params.Ki_q;

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

    v_int_d = 0;
    v_int_q = 0;

    t_hfi = 0;

end


% =========================================================================
% Protección básica
% =========================================================================

Imax = max(Imax, eps);
Ts   = max(Ts, eps);


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
% 2. Referencia de corriente q con compensación
% =========================================================================

% La compensación se interpreta directamente
% como corriente equivalente.

Iq_ref_raw = Iq_ref_req - Iq_comp;


% =========================================================================
% 3. Saturación vectorial de referencias
% =========================================================================

I_ref_mag_raw = hypot(Id_ref, Iq_ref_raw);

if I_ref_mag_raw > Imax

    scale_I = Imax/I_ref_mag_raw;

    Id_ref_lim = Id_ref     * scale_I;
    Iq_ref_lim = Iq_ref_raw * scale_I;

else

    Id_ref_lim = Id_ref;
    Iq_ref_lim = Iq_ref_raw;

end


% =========================================================================
% 4. Errores de corriente
% =========================================================================

err_d = Id_ref_lim - Id;
err_q = Iq_ref_lim - Iq;

error_d = err_d;
error_q = err_q;


% =========================================================================
% 5. Velocidad eléctrica
% =========================================================================

We = P * Wm;


% =========================================================================
% 6. Desacople y feedforward
%
% Modelo:
%
% Vd = R*Id + Ld*dId/dt - We*Lq*Iq
%
% Vq = R*Iq + Lq*dIq/dt
%      + We*Ld*Id
%      + We*Psi
%
% Se utilizan las referencias limitadas.
% =========================================================================

V_cross_d = ...
    -We * Lq * Iq_ref_lim;

V_cross_q = ...
     We * Ld * Id_ref_lim ...
   + We * Psi;


% =========================================================================
% 7. PI sin saturar
% =========================================================================

Vd_unsat = ...
    Kp_d*err_d ...
    + v_int_d ...
    + V_cross_d;

Vq_unsat = ...
    Kp_q*err_q ...
    + v_int_q ...
    + V_cross_q;


% =========================================================================
% 8. Límite de tensión disponible
% =========================================================================

% Actualmente se considera Vdc como máximo vectorial.
%
% Dependiendo del inversor/modulación podría cambiarse posteriormente a:
%
% Vdc/sqrt(3)
%
% si corresponde.

Vmax_inv = Vdc;


% Reserva de tensión para HFI
if HFI_enable == 1

    Vmax_control = ...
        max(Vmax_inv - abs(Amplitud_HFI), 0);

else

    Vmax_control = Vmax_inv;

end


% =========================================================================
% 9. Saturación vectorial de tensión fundamental
% =========================================================================

V_mag_unsat = hypot(Vd_unsat, Vq_unsat);

if V_mag_unsat > Vmax_control && ...
        V_mag_unsat > eps

    scale_V = ...
        Vmax_control/V_mag_unsat;

    Vd = Vd_unsat * scale_V;
    Vq = Vq_unsat * scale_V;

else

    Vd = Vd_unsat;
    Vq = Vq_unsat;

end


% =========================================================================
% 10. Anti-windup por back-calculation
% =========================================================================

Kaw_d = ...
    Ki_d/max(Kp_d, eps);

Kaw_q = ...
    Ki_q/max(Kp_q, eps);


v_int_d = ...
    v_int_d ...
    + Ts*( ...
        Ki_d*err_d ...
        + Kaw_d*(Vd - Vd_unsat));


v_int_q = ...
    v_int_q ...
    + Ts*( ...
        Ki_q*err_q ...
        + Kaw_q*(Vq - Vq_unsat));


% =========================================================================
% 11. Inyección HFI
% =========================================================================

t_hfi = t_hfi + Ts;


% Evita crecimiento ilimitado de t_hfi
if wh > 0

    T_hfi = 2*pi/wh;

    if t_hfi >= T_hfi

        t_hfi = t_hfi - T_hfi;

    end

end


if HFI_enable == 1

    v_hfi = ...
        Amplitud_HFI * sin(wh*t_hfi);

else

    v_hfi = 0;

end


% Inyección sobre el eje d estimado

Vd_hfi = Vd + v_hfi;
Vq_hfi = Vq;


% =========================================================================
% 12. Saturación final incluyendo HFI
% =========================================================================

V_mag_hfi = hypot(Vd_hfi, Vq_hfi);

if V_mag_hfi > Vmax_inv && ...
        V_mag_hfi > eps

    scale_hfi = ...
        Vmax_inv/V_mag_hfi;

    Vd_hfi = Vd_hfi * scale_hfi;
    Vq_hfi = Vq_hfi * scale_hfi;

end


% =========================================================================
% 13. Transformada inversa de Park
% =========================================================================

Valpha = ...
    Vd_hfi*cos_th ...
    - Vq_hfi*sin_th;

Vbetha = ...
    Vd_hfi*sin_th ...
    + Vq_hfi*cos_th;


% =========================================================================
% 14. Salidas de monitoreo
% =========================================================================

Id_out = Id;
Iq_out = Iq;

end