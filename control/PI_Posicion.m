function [Wm_ref, theta_ref_traj, e_theta, mode] = ...
    PI_Posicion(theta_target, theta_m, parameters)

% =========================================================================
% GENERADOR DE TRAYECTORIA PARA CONTROL DE POSICIÓN
%
% Entradas:
%   theta_target    : posición mecánica final deseada [rad]
%   theta_m         : posición mecánica medida/estimada [rad]
%   parameters      : estructura de parámetros
%
% Parámetros utilizados:
%   parameters.Ts
%   parameters.Wmax_pos
%   parameters.Amax_pos
%   parameters.theta_tol
%
% Salidas:
%   Wm_ref          : referencia de velocidad [rad/s]
%   theta_ref_traj  : referencia de posición de trayectoria [rad]
%   e_theta         : error de posición [rad]
%   mode            : 0 detenido
%                     1 acelerando
%                     2 velocidad constante
%                     3 frenando
%
% =========================================================================


% =========================================================================
% 0. Parámetros
% =========================================================================

Ts        = parameters.Ts;
Wmax      = parameters.Wmax_pos;
Amax      = parameters.Amax_pos;
theta_tol = parameters.theta_tol;


% Protección básica
Ts   = max(Ts, eps);
Wmax = max(abs(Wmax), eps);
Amax = max(abs(Amax), eps);
theta_tol = max(abs(theta_tol), 0);


% =========================================================================
% Estados persistentes
% =========================================================================

persistent Wm_cmd theta_traj

if isempty(Wm_cmd)

    Wm_cmd     = 0;
    theta_traj = theta_m;

end


% =========================================================================
% 1. Error de posición
% =========================================================================

e_theta = ...
    theta_target - theta_m;

dist = abs(e_theta);


% =========================================================================
% 2. Posición alcanzada
% =========================================================================

if dist <= theta_tol

    Wm_cmd = 0;

    theta_traj = theta_target;

    Wm_ref = 0;

    theta_ref_traj = theta_traj;

    mode = 0;

    return

end


% =========================================================================
% 3. Dirección del movimiento
% =========================================================================

dir = sign(e_theta);

if dir == 0
    dir = 1;
end


% =========================================================================
% 4. Distancia necesaria para frenar
%
%           W^2
% d = -------------
%         2*Amax
%
% =========================================================================

theta_brake = ...
    (Wm_cmd^2)/(2*Amax);


% =========================================================================
% 5. Perfil de velocidad
% =========================================================================

if dist <= theta_brake

    % ---------------------------------------------------------------------
    % Frenado
    % ---------------------------------------------------------------------

    mode = 3;

    Wmag = abs(Wm_cmd);

    Wmag = ...
        max( ...
            0, ...
            Wmag - Amax*Ts);

    Wm_cmd = ...
        dir * Wmag;


else

    % ---------------------------------------------------------------------
    % Aceleración / velocidad constante
    % ---------------------------------------------------------------------

    Wmag = abs(Wm_cmd);

    if Wmag < Wmax

        % Aceleración
        mode = 1;

        Wmag = ...
            min( ...
                Wmax, ...
                Wmag + Amax*Ts);

    else

        % Velocidad constante
        mode = 2;

        Wmag = Wmax;

    end


    Wm_cmd = ...
        dir * Wmag;

end


% =========================================================================
% 6. Integración de la trayectoria de posición
% =========================================================================

theta_traj = ...
    theta_traj ...
    + Wm_cmd*Ts;


% =========================================================================
% 7. Evitar sobrepasar la posición objetivo de trayectoria
% =========================================================================

if dir > 0

    if theta_traj > theta_target
        theta_traj = theta_target;
    end

else

    if theta_traj < theta_target
        theta_traj = theta_target;
    end

end


% =========================================================================
% 8. Salidas
% =========================================================================

Wm_ref = Wm_cmd;

theta_ref_traj = theta_traj;


end