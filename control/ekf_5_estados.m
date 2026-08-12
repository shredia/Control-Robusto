function x_out = ekf_5_estados(U, X, params)
% =========================================================================
% EKF 5 ESTADOS TIGHTLY COUPLED CON HFI
%
% Estados:
%
%   x = [Id; Iq; Wm; theta_e; Tx]
%
% Entradas:
%
%   U = [Va; Vb]
%   X = [Ia; Ib]
%
% Parámetros:
%
%   params.Ts
%   params.Ke
%   params.Kt
%   params.R
%   params.J
%   params.Ld
%   params.Lq
%   params.PP
%   params.B
%
%   params.wh
%
%   params.b0
%   params.b1
%   params.b2
%   params.a1
%   params.a2
%
%   params.HFI_flag
%
% Salida:
%
%   x_out = [Id; Iq; Wm; theta_e; Tx]
%
% =========================================================================


% =========================================================================
% 0. LECTURA DE PARÁMETROS
% =========================================================================

Ts = params.Ts;

Ke = params.Ke;
Kt = params.Kt;

R  = params.R;
J  = params.J;

Ld = params.Ld;
Lq = params.Lq;

PP = params.P;
B  = params.B;

wh = params.wh;

b0 = params.b0_hfi;
b1 = params.b1_hfi;
b2 = params.b2_hfi;

a1 = params.a1_hfi;
a2 = params.a2_hfi;

HFI_flag = params.HFI_flag;
Id_max = params.Id_max;
Iq_max = params.Iq_max;
Wm_max = params.Wm_max;
Tx_max = params.Tx_max;


% Protección numérica
Ts_eff = max(Ts, 1e-9);


% =========================================================================
% 1. ESTADOS PERSISTENTES
% =========================================================================

persistent x_hat
persistent Qk
persistent P
persistent xa1 xa2 ya1 ya2
persistent xb1 xb2 yb1 yb2

persistent t_h
persistent eps_f


if isempty(x_hat)

    % ---------------------------------------------------------------------
    % Estado inicial
    % ---------------------------------------------------------------------

    x_hat = zeros(5,1);


    % ---------------------------------------------------------------------
    % Covarianza inicial
    % ---------------------------------------------------------------------

    P = diag([ ...
        1e-6, ...   % Id
        1e-6, ...   % Iq
        1e-6, ...   % Wm
        1e-9, ...   % theta_e
        1e-5]);     % Tx


    % ---------------------------------------------------------------------
    % Ruido del proceso
    % ---------------------------------------------------------------------

    Qk = diag([ ...
        1e-8, ...    % Id
        1e-8, ...    % Iq
        1e-8, ...    % Wm
        1e-12, ...   % theta_e
        1e-5]);      % Tx


    % ---------------------------------------------------------------------
    % Estados del BPF de Ia
    % ---------------------------------------------------------------------

    xa1 = 0;
    xa2 = 0;

    ya1 = 0;
    ya2 = 0;


    % ---------------------------------------------------------------------
    % Estados del BPF de Ib
    % ---------------------------------------------------------------------

    xb1 = 0;
    xb2 = 0;

    yb1 = 0;
    yb2 = 0;


    % ---------------------------------------------------------------------
    % Estados HFI
    % ---------------------------------------------------------------------

    t_h   = 0;
    eps_f = 0;

end


% =========================================================================
% 2. MEDICIONES Y ENTRADAS
% =========================================================================

Ia_raw = X(1);
Ib_raw = X(2);

Va = U(1);
Vb = U(2);


% =========================================================================
% 3. ESTADOS ACTUALES
% =========================================================================

Id = x_hat(1);
Iq = x_hat(2);

Wm = x_hat(3);

% El estado angular se conserva continuo internamente
theta_e_state = x_hat(4);

% Torque desconocido
Tx = x_hat(5);


% Ángulo envuelto solamente para trigonometría
theta_e = wrap_pi(theta_e_state);

c = cos(theta_e);
s = sin(theta_e);


% =========================================================================
% 4. TRANSFORMACIÓN DE TENSIONES AL MARCO dq
% =========================================================================

Vd =  Va*c + Vb*s;
Vq = -Va*s + Vb*c;


% =========================================================================
% 5. PREDICCIÓN DEL ESTADO
% =========================================================================

We = PP*Wm;


% -------------------------------------------------------------------------
% Corriente d
% -------------------------------------------------------------------------

dId = ( ...
      Vd ...
    - R*Id ...
    + We*Lq*Iq) / Ld;


% -------------------------------------------------------------------------
% Corriente q
% -------------------------------------------------------------------------

dIq = ( ...
      Vq ...
    - R*Iq ...
    - We*Ld*Id ...
    - Ke*Wm) / Lq;


% -------------------------------------------------------------------------
% Torque electromagnético
% -------------------------------------------------------------------------
Krel = PP * (Lq-Ld);
Te = ...
      Kt*Iq ...
    - Krel*Id*Iq;


% -------------------------------------------------------------------------
% Dinámica mecánica
%
% Tx representa la carga desconocida
% -------------------------------------------------------------------------

dWm = ...
    (Te - B*Wm - Tx)/J;


% -------------------------------------------------------------------------
% Predicción
% -------------------------------------------------------------------------

x_pred = zeros(5,1);

x_pred(1) = Id + Ts_eff*dId;

x_pred(2) = Iq + Ts_eff*dIq;

x_pred(3) = Wm + Ts_eff*dWm;

% Ángulo eléctrico continuo
x_pred(4) = ...
    theta_e_state ...
    + Ts_eff*PP*Wm;

% Torque desconocido como random walk
x_pred(5) = Tx;


% =========================================================================
% 6. JACOBIANA DE TRANSICIÓN
% =========================================================================


A = [ ...

    -R/Ld, ...
     PP*Wm*Lq/Ld, ...
     PP*Lq*Iq/Ld, ...
     Vq/Ld, ...
     0;

    -PP*Wm*Ld/Lq, ...
    -R/Lq, ...
    -(PP*Ld*Id + Ke)/Lq, ...
    -Vd/Lq, ...
     0;

     -Krel*Iq/J, ...
     (Kt - Krel*Id)/J, ...
    -B/J, ...
     0, ...
    -1/J;

     0, ...
     0, ...
     PP, ...
     0, ...
     0;

     0, ...
     0, ...
     0, ...
     0, ...
     0];


Phi = eye(5) + Ts_eff*A;


% =========================================================================
% 7. PREDICCIÓN DE COVARIANZA
% =========================================================================

P_pred = ...
    Phi*P*Phi' ...
    + Qk;

% Forzar simetría numérica
P_pred = ...
    0.5*(P_pred + P_pred');


% =========================================================================
% 8. FILTRO PASABANDA DE LAS CORRIENTES HFI
% =========================================================================

Iah = ...
      b0*Ia_raw ...
    + b1*xa1 ...
    + b2*xa2 ...
    - a1*ya1 ...
    - a2*ya2;


xa2 = xa1;
xa1 = Ia_raw;

ya2 = ya1;
ya1 = Iah;


Ibh = ...
      b0*Ib_raw ...
    + b1*xb1 ...
    + b2*xb2 ...
    - a1*yb1 ...
    - a2*yb2;


xb2 = xb1;
xb1 = Ib_raw;

yb2 = yb1;
yb1 = Ibh;


% =========================================================================
% 9. PROYECCIÓN HFI RESPECTO DEL ÁNGULO PREDICHO
% =========================================================================

theta_e_pred = wrap_pi(x_pred(4));

ch = cos(theta_e_pred);
sh = sin(theta_e_pred);


% Componente HFI sobre el eje q estimado
Iqh = ...
    -Iah*sh ...
    + Ibh*ch;


% =========================================================================
% 10. DEMODULACIÓN SÍNCRONA
% =========================================================================

if wh > 0

    t_h = t_h + Ts_eff;

    T_hfi = 2*pi/wh;

    if t_h >= T_hfi

        t_h = ...
            t_h ...
            - T_hfi*floor(t_h/T_hfi);

    end


    % IMPORTANTE:
    %
    % Verificar signo y fase según la señal
    % efectivamente inyectada.
    %
    eps_raw = ...
        -Iqh*cos(wh*t_h);

else

    t_h = 0;

    eps_raw = 0;

end


% =========================================================================
% 11. FILTRO PASABAJOS DEL ERROR HFI
% =========================================================================

f_lpf_hfi = 50;

tau_lpf = ...
    1/(2*pi*f_lpf_hfi);

alpha_lpf = ...
    Ts_eff/(Ts_eff + tau_lpf);


eps_f = ...
    eps_f ...
    + alpha_lpf*(eps_raw - eps_f);


% =========================================================================
% 12. NORMALIZACIÓN DEL ERROR HFI
% =========================================================================

DeltaL = ...
    Lq - Ld;

Lavg = ...
    0.5*(Ld + Lq);


if abs(DeltaL) < 1e-12 || Lavg <= 0

    eps_norm = 0;

else

    saliency = ...
        DeltaL/Lavg;


    % ---------------------------------------------------------------------
    % Normalización original
    % ---------------------------------------------------------------------

    eps_norm = ...
        eps_f/abs(saliency);


    % ---------------------------------------------------------------------
    % Límite de seguridad
    %
    % ±45 grados eléctricos
    % ---------------------------------------------------------------------

    max_eps = pi/4;

    eps_norm = ...
        max( ...
            -max_eps, ...
            min(max_eps, eps_norm));

end


% -------------------------------------------------------------------------
% Desactivar HFI
% -------------------------------------------------------------------------

if ~HFI_flag

    eps_norm = 0;

end


% -------------------------------------------------------------------------
% Protección numérica
% -------------------------------------------------------------------------

if ~isfinite(eps_norm)

    eps_norm = 0;

end


% =========================================================================
% 13. PREDICCIÓN DE LAS MEDICIONES DE CORRIENTE
% =========================================================================

Id_p = x_pred(1);
Iq_p = x_pred(2);

cp = cos(theta_e_pred);
sp = sin(theta_e_pred);


Ia_pred = ...
      Id_p*cp ...
    - Iq_p*sp;


Ib_pred = ...
      Id_p*sp ...
    + Iq_p*cp;


% =========================================================================
% 14. VECTOR DE MEDICIÓN
% =========================================================================
%
% La tercera medición impone:
%
% eps_norm -> 0
%
% =========================================================================

y_med = [ ...
    Ia_raw;
    Ib_raw;
    0];


y_pred = [ ...
    Ia_pred;
    Ib_pred;
    eps_norm];


innov = ...
    y_med - y_pred;


% =========================================================================
% 15. JACOBIANA DE MEDICIÓN
% =========================================================================

dIa_dtheta_e = ...
    -Id_p*sp ...
    -Iq_p*cp;


dIb_dtheta_e = ...
     Id_p*cp ...
    -Iq_p*sp;


% Si:
%
% eps_norm ~= theta_real - theta_estimado
%
% entonces:
%
% d(eps_norm)/d(theta_estimado) ~= -1
%
dEps_dtheta_e = -1;


Ck = [ ...

    cp, -sp, 0, dIa_dtheta_e, 0;

    sp,  cp, 0, dIb_dtheta_e, 0;

    0,   0,  0, dEps_dtheta_e, 0];


% =========================================================================
% 16. COVARIANZA DE LAS MEDICIONES
% =========================================================================

sigma_Ia = 0.1;
sigma_Ib = 0.1;


if HFI_flag

    sigma_HFI = 0.25;

else

    % HFI prácticamente ignorado
    sigma_HFI = 1e6;

end


Rk = diag([ ...
    sigma_Ia^2, ...
    sigma_Ib^2, ...
    sigma_HFI^2]);


% =========================================================================
% 17. GANANCIA DE KALMAN
% =========================================================================

Sk = ...
    Ck*P_pred*Ck' ...
    + Rk;


% Forzar simetría
Sk = ...
    0.5*(Sk + Sk');


% Evitar inversión explícita
Lk = ...
    (P_pred*Ck')/Sk;


% =========================================================================
% 18. CORRECCIÓN TIGHTLY COUPLED
% =========================================================================

x_corr_ekf = ...
    x_pred ...
    + Lk*innov;


% El ángulo NO se envuelve internamente.
%
% x_corr_ekf(4) permanece continuo.


% =========================================================================
% 19. ACTUALIZACIÓN DE COVARIANZA
%
% Forma de Joseph
% =========================================================================

I5 = eye(5);


P = ...
    (I5 - Lk*Ck) ...
    *P_pred ...
    *(I5 - Lk*Ck)' ...
    + Lk*Rk*Lk';


% Forzar simetría
P = ...
    0.5*(P + P');


% =========================================================================
% 20. PROTECCIÓN NUMÉRICA
% =========================================================================

if any(~isfinite(x_corr_ekf)) || ...
   any(~isfinite(P(:)))

    x_corr_ekf = x_pred;

    P = P_pred;

end


% =========================================================================
% 21. ACTUALIZACIÓN DEL ESTADO
% =========================================================================

x_hat = x_corr_ekf;
x_hat(1) = min(max(x_hat(1), -Id_max), Id_max);
x_hat(2) = min(max(x_hat(2), -Iq_max), Iq_max);
x_hat(3) = min(max(x_hat(3), -Wm_max), Wm_max);
x_hat(5) = min(max(x_hat(5), -Tx_max), Tx_max);
% =========================================================================
% 22. SALIDA
% =========================================================================

x_out = x_hat;


% Solo se envuelve el ángulo de salida
x_out(4) = wrap_pi(x_hat(4));


end


% =========================================================================
% FUNCIÓN AUXILIAR
% =========================================================================

function out = wrap_pi(in)

out = atan2(sin(in), cos(in));

end