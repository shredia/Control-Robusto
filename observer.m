function observer(block)
% Level-2 MATLAB S-Function (continuous-time)
% Observador Luenberger no lineal en alpha-beta
% Inports:
%   1) u  = [vA; vB]
%   2) y  = [iA_meas; iB_meas]
%   3) Lc = ganancia del observador (5x2) o vector 10x1
% Outport:
%   xhat = [iA; iB; we; th_wrapped; TL]

  setup(block);
end

%==================== CONFIG (edita si quieres embebido) ====================%
function [P, x0] = get_config()
  % Parámetros físicos
  P.R  = 3.4;      % Ohm
  P.L  = 6e-3;     % H
  P.Ke = 5.5e-3;   % V*s/rad
  P.p  = 50;       % pares de polos
  P.Kt = P.Ke*P.p; % N*m/A (convención usada aquí)
  P.B  = 5e-5;     % N*m*s/rad
  P.J  = 8e-3;     % kg*m^2

  % Estado inicial
  x0 = [0;0;0;0;0];
end
%============================================================================%

function setup(block)
  % Puertos
  block.NumInputPorts  = 3;
  block.NumOutputPorts = 1;

  % Inport 1: u = [vA; vB]
  block.SetPreCompInpPortInfoToDynamic;
  block.InputPort(1).Dimensions        = 2;
  block.InputPort(1).DatatypeID        = 0;
  block.InputPort(1).Complexity        = 'Real';
  block.InputPort(1).DirectFeedthrough = true;

  % Inport 2: y = [iA_meas; iB_meas]
  block.InputPort(2).Dimensions        = 2;
  block.InputPort(2).DatatypeID        = 0;
  block.InputPort(2).Complexity        = 'Real';
  block.InputPort(2).DirectFeedthrough = true;

% Inport 3: Lc (matriz 5x2)
block.InputPort(3).Dimensions        = [5 2];
block.InputPort(3).DatatypeID        = 0;
block.InputPort(3).Complexity        = 'Real';
block.InputPort(3).DirectFeedthrough = true;

  % Outport: xhat (5x1)
  block.SetPreCompOutPortInfoToDynamic;
  block.OutputPort(1).Dimensions       = 5;
  block.OutputPort(1).DatatypeID       = 0;
  block.OutputPort(1).Complexity       = 'Real';

  % Estados continuos y tiempo continuo
  block.NumContStates = 5;
  block.SampleTimes   = [0 0];

  % ¡No declares NumDworks aquí!
  % Métodos
  block.SimStateCompliance = 'DefaultSimState';
  block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions', @InitConditions);
  block.RegBlockMethod('Outputs',              @Outputs);
  block.RegBlockMethod('Derivatives',          @Derivatives);
end

function DoPostPropSetup(block)
  block.NumDworks = 1;

  block.Dwork(1).Name            = 'Pvec';
  block.Dwork(1).Dimensions      = 7;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;
end

function InitConditions(block)
  [P, x0] = get_config();             % ✅ cargar aquí
  block.Dwork(1).Data = [P.R, P.L, P.Ke, P.Kt, P.B, P.J, P.p];
  block.ContStates.Data = x0(:);
end

function Outputs(block)
  x = block.ContStates.Data; % [iA;iB;we;th;TL]
  x_out = x;
  x_out(4) = atan2(sin(x(4)), cos(x(4))); % wrap solo para mostrar
  block.OutputPort(1).Data = x_out;
end

function Derivatives(block)
  % Parámetros físicos
  Pvec = block.Dwork(1).Data(:).';
  Rm = Pvec(1); Lm = Pvec(2); Ke = Pvec(3);
  Kt = Pvec(4); Bm = Pvec(5); J  = Pvec(6); p = Pvec(7);
  gamma = p / J;

  % Entradas
  u  = block.InputPort(1).Data; % [vA; vB]
  y  = block.InputPort(2).Data; % [iA_meas; iB_meas]
  Lc = block.InputPort(3).Data;

  % Validación opcional
if ~isequal(size(Lc), [5 2])
    error('El puerto 3 debe ser una matriz 5x2.');
end

  % Estados
  x  = block.ContStates.Data;   % [iA;iB;we;th;TL]
  iA = x(1); iB = x(2); we = x(3); th = x(4); TL = x(5);
  vA = u(1); vB = u(2);

  % Dinámica no lineal f(x,u) en alpha-beta
  diA = (vA - Rm*iA + we*Ke*sin(th)) / Lm;
  diB = (vB - Rm*iB - we*Ke*cos(th)) / Lm;
  dwe = gamma*( Kt*(cos(th)*iB - sin(th)*iA) - TL ) - gamma*Bm*we;
  dth = we;
  dTL = 0;

  fxu = [diA; diB; dwe; dth; dTL];

  % Residuo (y - C*xhat), C=[1 0 0 0 0; 0 1 0 0 0]
  res = [ y(1) - iA;
          y(2) - iB ];

  % Observador continuo
  xdot = fxu + Lc * res;

  block.Derivatives.Data = xdot;
end
