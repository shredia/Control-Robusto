function [T_res, error_res] = Resonant_Fixed(Wm_ref, Wm_filt, params)
% =========================================================================
% CONTROL RESONANTE DE VELOCIDAD
%
% Genera una referencia adicional de torque para compensar una
% perturbación periódica, por ejemplo cogging.
%
% Después:
%
%       T_ref = T_PI + T_res
%
% y T_ref se entrega al algoritmo MTPA.
%
% =========================================================================

persistent x1 x2

if isempty(x1)
    x1 = 0.0;
    x2 = 0.0;
end

Ts       = params.Ts;
Kr       = params.Kr;
zeta     = params.zeta;
wr       = params.wr;
T_res_max = params.T_res_max;
enable   = params.enable;

% Error de velocidad
error_res = Wm_ref - Wm_filt;

% -------------------------------------------------------------------------
% Desactivar resonador
% -------------------------------------------------------------------------

if enable == 0

    x1 = 0.0;
    x2 = 0.0;

    T_res = 0.0;

    return;
end

% -------------------------------------------------------------------------
% Resonador
%
%                 2*zeta*wr*s
% G(s) = Kr ---------------------------
%             s² + 2*zeta*wr*s + wr²
%
% -------------------------------------------------------------------------

dx1 = x2;

dx2 = ...
    -(wr * wr) * x1 ...
    -2.0 * zeta * wr * x2 ...
    +2.0 * zeta * wr * Kr * error_res;

% Euler discreto
x1 = x1 + Ts * dx1;
x2 = x2 + Ts * dx2;

% Salida en torque
T_res = x2;

% Saturación
if T_res > T_res_max
    T_res = T_res_max;

elseif T_res < -T_res_max
    T_res = -T_res_max;
end

end