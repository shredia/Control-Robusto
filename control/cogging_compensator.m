function [Iq_cog, e_h, y1_out, y2_out, y3_out, ...
          f1_out, f2_out, f3_out] = ...
    cogging_compensator( ...
    Wm_hat, Wm_filt, Ts, ...
    Nc, ...
    K1, K2, K3, ...
    Q1, Q2, Q3, ...
    Iq_cog_max, ...
    Wm_min, ...
    enable)

% =========================================================================
% COMPENSADOR ARMÓNICO ADAPTATIVO DE COGGING
%
% Extracción de 3 armónicos del ripple de velocidad:
%
%   e_h = Wm_hat - Wm_filt
%
% Frecuencia fundamental de cogging:
%
%   f1 = Nc * abs(Wm_filt)/(2*pi)
%
% Armónicos:
%
%   f2 = 2*f1
%   f3 = 3*f1
%
% Compensación:
%
%   Iq_cog = K1*y1 + K2*y2 + K3*y3
%
%
% ENTRADAS
% -------------------------------------------------------------------------
% Wm_hat       : velocidad estimada EKF [rad/s]
% Wm_filt      : velocidad filtrada [rad/s]
% Ts           : tiempo de muestreo [s]
%
% Nc           : periodicidad fundamental por vuelta mecánica
%
% K1,K2,K3     : ganancias de compensación [A/(rad/s)]
%
% Q1,Q2,Q3     : factor de calidad de cada BPF
%
% Iq_cog_max   : saturación TOTAL de compensación [A]
%
% Wm_min       : velocidad mínima para activar compensación [rad/s]
%
% enable       : 1 = activo
%                0 = desactivado
%
%
% SALIDAS
% -------------------------------------------------------------------------
% Iq_cog       : corriente total de compensación [A]
%
% e_h          : ripple de velocidad [rad/s]
%
% y1_out       : componente armónica fundamental [rad/s]
% y2_out       : segundo armónico [rad/s]
% y3_out       : tercer armónico [rad/s]
%
% f1_out       : frecuencia fundamental [Hz]
% f2_out       : segundo armónico [Hz]
% f3_out       : tercer armónico [Hz]
%
% =========================================================================


% =========================================================================
% ESTADOS PERSISTENTES
%
% Cada biquad necesita:
%
% x[k-1]
% x[k-2]
% y[k-1]
% y[k-2]
% =========================================================================

persistent x11 x12 yy11 yy12
persistent x21 x22 yy21 yy22
persistent x31 x32 yy31 yy32

if isempty(x11)

    x11  = 0.0;
    x12  = 0.0;
    yy11 = 0.0;
    yy12 = 0.0;

    x21  = 0.0;
    x22  = 0.0;
    yy21 = 0.0;
    yy22 = 0.0;

    x31  = 0.0;
    x32  = 0.0;
    yy31 = 0.0;
    yy32 = 0.0;

end


% =========================================================================
% INICIALIZACIÓN EXPLÍCITA DE SALIDAS
% =========================================================================

Iq_cog = 0.0;

e_h = 0.0;

y1_out = 0.0;
y2_out = 0.0;
y3_out = 0.0;

f1_out = 0.0;
f2_out = 0.0;
f3_out = 0.0;


% =========================================================================
% PROTECCIONES
% =========================================================================

Ts_eff = max(Ts,1e-9);

fs = 1/Ts_eff;

Nc_eff = abs(Nc);

Q1_eff = max(Q1,0.5);
Q2_eff = max(Q2,0.5);
Q3_eff = max(Q3,0.5);

Iq_max_eff = abs(Iq_cog_max);

Wm_min_eff = abs(Wm_min);


% =========================================================================
% 1. COMPONENTE RÁPIDA DE VELOCIDAD
% =========================================================================

e_h = Wm_hat - Wm_filt;


% =========================================================================
% 2. FRECUENCIA FUNDAMENTAL DEL COGGING
%
% IMPORTANTE:
%
% Se usa Wm_filt para calcular la frecuencia.
%
% Esto evita que el propio ripple module continuamente
% la frecuencia central de los filtros.
% =========================================================================

f_mech = abs(Wm_filt)/(2*pi);

f1 = Nc_eff*f_mech;

f2 = 2.0*f1;

f3 = 3.0*f1;


% =========================================================================
% 3. LÍMITE DE FRECUENCIA
%
% Dejamos margen respecto de Nyquist.
% =========================================================================

f_max = 0.40*fs;


% =========================================================================
% 4. DESACTIVAR A MUY BAJA VELOCIDAD
% =========================================================================

if (enable < 0.5) || (abs(Wm_filt) < Wm_min_eff)

    % -------------------------------------------------------------
    % Reset suave de estados
    % -------------------------------------------------------------

    decay = 0.95;

    x11  = decay*x11;
    x12  = decay*x12;
    yy11 = decay*yy11;
    yy12 = decay*yy12;

    x21  = decay*x21;
    x22  = decay*x22;
    yy21 = decay*yy21;
    yy22 = decay*yy22;

    x31  = decay*x31;
    x32  = decay*x32;
    yy31 = decay*yy31;
    yy32 = decay*yy32;

    Iq_cog = 0.0;

    y1_out = 0.0;
    y2_out = 0.0;
    y3_out = 0.0;

    f1_out = 0.0;
    f2_out = 0.0;
    f3_out = 0.0;

    return;

end


% =========================================================================
% 5. LIMITAR FRECUENCIAS INDIVIDUALES
% =========================================================================

if f1 > f_max
    f1 = f_max;
end

if f2 > f_max
    f2 = f_max;
end

if f3 > f_max
    f3 = f_max;
end


f1_out = f1;
f2_out = f2;
f3_out = f3;


% =========================================================================
% 6. ARMÓNICO 1
% =========================================================================

w1 = 2*pi*f1/fs;

alpha1 = sin(w1)/(2*Q1_eff);

a01 = 1 + alpha1;

b01 =  alpha1/a01;
b11 =  0.0;
b21 = -alpha1/a01;

a11 = -2*cos(w1)/a01;
a21 = (1-alpha1)/a01;


y1 = ...
      b01*e_h ...
    + b11*x11 ...
    + b21*x12 ...
    - a11*yy11 ...
    - a21*yy12;


% Actualizar estados

x12 = x11;
x11 = e_h;

yy12 = yy11;
yy11 = y1;


% =========================================================================
% 7. ARMÓNICO 2
% =========================================================================

w2 = 2*pi*f2/fs;

alpha2 = sin(w2)/(2*Q2_eff);

a02 = 1 + alpha2;

b02 =  alpha2/a02;
b12 =  0.0;
b22 = -alpha2/a02;

a12 = -2*cos(w2)/a02;
a22 = (1-alpha2)/a02;


y2 = ...
      b02*e_h ...
    + b12*x21 ...
    + b22*x22 ...
    - a12*yy21 ...
    - a22*yy22;


% Actualizar estados

x22 = x21;
x21 = e_h;

yy22 = yy21;
yy21 = y2;


% =========================================================================
% 8. ARMÓNICO 3
% =========================================================================

w3 = 2*pi*f3/fs;

alpha3 = sin(w3)/(2*Q3_eff);

a03 = 1 + alpha3;

b03 =  alpha3/a03;
b13 =  0.0;
b23 = -alpha3/a03;

a13 = -2*cos(w3)/a03;
a23 = (1-alpha3)/a03;


y3 = ...
      b03*e_h ...
    + b13*x31 ...
    + b23*x32 ...
    - a13*yy31 ...
    - a23*yy32;


% Actualizar estados

x32 = x31;
x31 = e_h;

yy32 = yy31;
yy31 = y3;


% =========================================================================
% 9. SALIDAS DE LOS ARMÓNICOS
% =========================================================================

y1_out = y1;
y2_out = y2;
y3_out = y3;


% =========================================================================
% 10. COMPENSACIÓN INDIVIDUAL
% =========================================================================

Iq1 = K1*y1;
Iq2 = K2*y2;
Iq3 = K3*y3;


% =========================================================================
% 11. SUMA DE ARMÓNICOS
% =========================================================================

Iq_unsat = Iq1 + Iq2 + Iq3;


% =========================================================================
% 12. SATURACIÓN TOTAL
% =========================================================================

if Iq_unsat > Iq_max_eff

    Iq_cog = Iq_max_eff;

elseif Iq_unsat < -Iq_max_eff

    Iq_cog = -Iq_max_eff;

else

    Iq_cog = Iq_unsat;

end


end