function resultado = FinishFcn(simOut, Wm_operacion, timeWindow, doPrint)

% =========================================================================
% FINISHFCN
%
% Analiza M1 y M2 a partir de un objeto Simulink.SimulationOutput.
%
% Entradas:
%
%   simOut
%       Simulink.SimulationOutput
%
%   Wm_operacion
%       Velocidad mecánica de operación [rad/s]
%
%   timeWindow
%       [t1 t2] ventana estacionaria
%
%   doPrint
%       true  -> imprime resumen
%       false -> no imprime
%
%
% Salida:
%
%   resultado.M1
%   resultado.M2
%
% Cada motor contiene:
%
%   .phi_Wm
%   .gain_Wm
%   .delay_Wm
%   .samples_Wm
%
%   .phi_Iq
%   .gain_Iq
%   .delay_Iq
%   .samples_Iq
%
%   .phi_theta
%   .gain_theta
%   .delay_theta
%   .samples_theta
%
%   .phi_TPR
%   .gain_TPR
%   .delay_TPR
%   .samples_TPR
%
%   .Wm_mean
%   .Wm_std
%   .Wm_ripple_pp
%   .Wm_ripple_amp
%   .Wm_RMSE
%   .Wm_error_mean
%
% resultado también contiene:
%
%   .Nc
%   .wcog
%   .fcog
%   .Wm_operacion
%   .timeWindow
%
% =========================================================================


%% ========================================================================
% VALORES POR DEFECTO
% ========================================================================

if nargin < 2 || isempty(Wm_operacion)

    Wm_operacion = 10;

end


if nargin < 3 || isempty(timeWindow)

    timeWindow = [0.12 0.20];

end


if nargin < 4 || isempty(doPrint)

    doPrint = true;

end


%% ========================================================================
% FRECUENCIA DE COGGING
% ========================================================================

Nc = 200;

Wm_operacion = abs(Wm_operacion);

wcog = ...
    Nc * Wm_operacion;

fcog = ...
    wcog / (2*pi);


resultado.Nc = Nc;

resultado.wcog = wcog;

resultado.fcog = fcog;

resultado.Wm_operacion = Wm_operacion;

resultado.timeWindow = timeWindow;


%% ========================================================================
% EXTRAER SEÑALES M1
% ========================================================================

M1.Wm_real = obtener_senal(simOut,'Wm_real_M1');
M1.Wm_obs  = obtener_senal(simOut,'Wm_obs_M1');

M1.Iq_ref = obtener_senal(simOut,'Iq_ref_M1');
M1.Iq     = obtener_senal(simOut,'Iq_M1');

M1.theta_real = obtener_senal(simOut,'theta_e_planta_M1');
M1.theta_obs  = obtener_senal(simOut,'theta_e_obs_M1');

M1.Tdtm = obtener_senal(simOut,'Tdtm_M1');

% Si T_res corresponde a la salida del resonador/PR
M1.T_PR = obtener_senal(simOut,'T_res_M1');


%% ========================================================================
% EXTRAER SEÑALES M2
% ========================================================================

M2.Wm_real = obtener_senal(simOut,'Wm_real_M2');
M2.Wm_obs  = obtener_senal(simOut,'Wm_obs_M2');

M2.Iq_ref = obtener_senal(simOut,'Iq_ref_M2');
M2.Iq     = obtener_senal(simOut,'Iq_M2');

M2.theta_real = obtener_senal(simOut,'theta_e_planta_M2');
M2.theta_obs  = obtener_senal(simOut,'theta_e_obs_M2');

M2.Tdtm = obtener_senal(simOut,'Tdtm_M2');

% Si T_res corresponde a la salida resonante
M2.T_PR = obtener_senal(simOut,'T_res_M2');

%% ========================================================================
% ANALIZAR MOTOR 1
% ========================================================================

resultado.M1 = ...
    analizar_motor( ...
        M1, ...
        fcog, ...
        Wm_operacion, ...
        timeWindow);


%% ========================================================================
% ANALIZAR MOTOR 2
% ========================================================================

resultado.M2 = ...
    analizar_motor( ...
        M2, ...
        fcog, ...
        Wm_operacion, ...
        timeWindow);


%% ========================================================================
% IMPRIMIR
% ========================================================================

if doPrint

    fprintf('\n');
    fprintf('=============================================================\n');
    fprintf(' ANÁLISIS DE SIMULACIÓN\n');
    fprintf('=============================================================\n');

    fprintf('Wm operación = %.4f rad/s\n',Wm_operacion);
    fprintf('wcog         = %.4f rad/s\n',wcog);
    fprintf('fcog         = %.4f Hz\n',fcog);

    fprintf( ...
        'Ventana      = %.4f - %.4f s\n', ...
        timeWindow(1), ...
        timeWindow(2));

    fprintf('=============================================================\n');


    imprimir_motor( ...
        resultado.M1, ...
        'M1');


    imprimir_motor( ...
        resultado.M2, ...
        'M2');


    fprintf('=============================================================\n');
    fprintf(' RESUMEN DE FASE\n');
    fprintf('=============================================================\n');

    fprintf( ...
        'Wm_obs/Wm_real        M1 : %+8.3f deg\n', ...
        resultado.M1.phi_Wm);

    fprintf( ...
        'Wm_obs/Wm_real        M2 : %+8.3f deg\n', ...
        resultado.M2.phi_Wm);

    fprintf( ...
        'Iq/Iq_ref             M1 : %+8.3f deg\n', ...
        resultado.M1.phi_Iq);

    fprintf( ...
        'Iq/Iq_ref             M2 : %+8.3f deg\n', ...
        resultado.M2.phi_Iq);

    fprintf( ...
        'theta_obs/theta_real  M1 : %+8.3f deg\n', ...
        resultado.M1.phi_theta);

    fprintf( ...
        'theta_obs/theta_real  M2 : %+8.3f deg\n', ...
        resultado.M2.phi_theta);

    fprintf( ...
        'T_PR/Tdtm             M1 : %+8.3f deg\n', ...
        resultado.M1.phi_TPR);

    fprintf( ...
        'T_PR/Tdtm             M2 : %+8.3f deg\n', ...
        resultado.M2.phi_TPR);

    fprintf('=============================================================\n');

end


end


%% ========================================================================
% FUNCIÓN: ANALIZAR UN MOTOR
% ========================================================================

function r = analizar_motor( ...
    S, ...
    fcog, ...
    Wm_operacion, ...
    timeWindow)


% =========================================================================
% Wm OBSERVADA / REAL
% =========================================================================

[ ...
    r.phi_Wm, ...
    r.gain_Wm, ...
    r.delay_Wm, ...
    r.samples_Wm ...
] = fase_timeseries( ...
        S.Wm_real, ...
        S.Wm_obs, ...
        fcog, ...
        timeWindow);


% =========================================================================
% Iq / Iq_ref
% =========================================================================

[ ...
    r.phi_Iq, ...
    r.gain_Iq, ...
    r.delay_Iq, ...
    r.samples_Iq ...
] = fase_timeseries( ...
        S.Iq_ref, ...
        S.Iq, ...
        fcog, ...
        timeWindow);


% =========================================================================
% theta OBSERVADA / REAL
% =========================================================================

[ ...
    r.phi_theta, ...
    r.gain_theta, ...
    r.delay_theta, ...
    r.samples_theta ...
] = fase_timeseries( ...
        S.theta_real, ...
        S.theta_obs, ...
        fcog, ...
        timeWindow);


% =========================================================================
% T_PR / Tdtm
% =========================================================================

[ ...
    r.phi_TPR, ...
    r.gain_TPR, ...
    r.delay_TPR, ...
    r.samples_TPR ...
] = fase_timeseries( ...
        S.Tdtm, ...
        S.T_PR, ...
        fcog, ...
        timeWindow);


% =========================================================================
% MÉTRICAS DE VELOCIDAD
% =========================================================================

t = ...
    S.Wm_real.Time(:);

Wm = ...
    squeeze(S.Wm_real.Data);

Wm = ...
    Wm(:);


idx = ...
    t >= timeWindow(1) & ...
    t <= timeWindow(2);


Wm = ...
    Wm(idx);


Wm = ...
    Wm(isfinite(Wm));


if isempty(Wm)

    error( ...
        'No existen datos de velocidad dentro de la ventana.');

end


r.Wm_mean = ...
    mean(Wm);


r.Wm_std = ...
    std(Wm);


r.Wm_ripple_pp = ...
    max(Wm) - min(Wm);


r.Wm_ripple_amp = ...
    r.Wm_ripple_pp/2;


r.Wm_RMSE = ...
    sqrt( ...
        mean( ...
            (Wm - Wm_operacion).^2));


r.Wm_error_mean = ...
    r.Wm_mean - Wm_operacion;


end


%% ========================================================================
% FUNCIÓN: OBTENER SEÑAL DESDE simOut
% ========================================================================
%% ========================================================================
% FUNCIÓN: OBTENER SEÑAL DESDE logsout
% ========================================================================

function x = obtener_senal(simOut, nombre)

% =========================================================================
% OBTENER_SENAL
%
% Busca exclusivamente una señal dentro del Dataset:
%
%       simOut.logsout
%
% Devuelve:
%
%       timeseries
%
% Esto permite mantener todas las señales utilizadas por FinishFcn
% agrupadas dentro de Signal Logging.
% =========================================================================


%% ========================================================================
% VERIFICAR SimulationOutput
% ========================================================================

if ~isa(simOut,'Simulink.SimulationOutput')

    error( ...
        'FinishFcn:InvalidSimulationOutput', ...
        'simOut debe ser un objeto Simulink.SimulationOutput.');

end


%% ========================================================================
% VERIFICAR logsout
% ========================================================================

variables = simOut.who;

if ~any(strcmp(variables,'logsout'))

    error( ...
        'FinishFcn:LogsoutNotFound', ...
        ['No existe "logsout" dentro de SimulationOutput. ' ...
         'Verifique que Signal Logging esté activado.']);

end


logsout = simOut.get('logsout');


if isempty(logsout)

    error( ...
        'FinishFcn:LogsoutEmpty', ...
        'El Dataset logsout está vacío.');

end


%% ========================================================================
% BUSCAR SEÑAL
% ========================================================================

try

    elemento = logsout.getElement(nombre);

catch

    elemento = [];

end


%% ========================================================================
% VERIFICAR RESULTADO
% ========================================================================

if ~isempty(elemento)

    x = elemento.Values;

    return;

end


%% ========================================================================
% SEÑAL NO ENCONTRADA
% ========================================================================

fprintf('\n');
fprintf('=============================================================\n');
fprintf(' SEÑAL NO ENCONTRADA EN logsout\n');
fprintf('=============================================================\n');

fprintf('Buscada: %s\n\n',nombre);

fprintf('Señales disponibles:\n');

for k = 1:logsout.numElements

    elem = logsout.getElement(k);

    fprintf( ...
        '%3d | %s\n', ...
        k, ...
        elem.Name);

end

fprintf('=============================================================\n');


error( ...
    'FinishFcn:SignalNotFound', ...
    'No se encontró la señal "%s" dentro de logsout.', ...
    nombre);

end
%% ========================================================================
% FUNCIÓN: IMPRIMIR RESULTADOS DE UN MOTOR
% ========================================================================

function imprimir_motor(r,nombre)


fprintf('\n');
fprintf('-------------------------------------------------------------\n');
fprintf(' MOTOR %s\n',nombre);
fprintf('-------------------------------------------------------------\n');


fprintf( ...
    'Wm_obs / Wm_real      : %+9.3f deg   Gain = %.4f\n', ...
    r.phi_Wm, ...
    r.gain_Wm);


fprintf( ...
    'Iq / Iq_ref           : %+9.3f deg   Gain = %.4f\n', ...
    r.phi_Iq, ...
    r.gain_Iq);


fprintf( ...
    'theta_obs/theta_real  : %+9.3f deg   Gain = %.4f\n', ...
    r.phi_theta, ...
    r.gain_theta);


fprintf( ...
    'T_PR / Tdtm           : %+9.3f deg   Gain = %.4f\n', ...
    r.phi_TPR, ...
    r.gain_TPR);


fprintf('-------------------------------------------------------------\n');


fprintf( ...
    'Wm media              : %.6f rad/s\n', ...
    r.Wm_mean);


fprintf( ...
    'Std Wm                : %.6f rad/s\n', ...
    r.Wm_std);


fprintf( ...
    'Ripple Wm p-p         : %.6f rad/s\n', ...
    r.Wm_ripple_pp);


fprintf( ...
    'Ripple Wm amplitud    : %.6f rad/s\n', ...
    r.Wm_ripple_amp);


fprintf( ...
    'RMSE Wm               : %.6f rad/s\n', ...
    r.Wm_RMSE);


fprintf( ...
    'Error medio Wm        : %+.6f rad/s\n', ...
    r.Wm_error_mean);


end


%% ========================================================================
% FUNCIÓN: FASE ENTRE DOS TIMESERIES
% ========================================================================

function [phi_deg, gain_ratio, delay_s, delay_samples] = ...
    fase_timeseries( ...
        x_ref, ...
        y, ...
        f0, ...
        timeWindow)


% =========================================================================
% EXTRAER DATOS
% =========================================================================

tx = ...
    x_ref.Time(:);

x = ...
    squeeze(x_ref.Data);

x = ...
    x(:);


ty = ...
    y.Time(:);

yd = ...
    squeeze(y.Data);

yd = ...
    yd(:);


% =========================================================================
% INTERPOLAR
% =========================================================================

if length(tx) ~= length(ty) || ...
        any(abs(tx-ty) > 1e-12)

    yd = ...
        interp1( ...
            ty, ...
            yd, ...
            tx, ...
            'linear', ...
            NaN);

end


% =========================================================================
% VENTANA
% =========================================================================

idx = ...
    tx >= timeWindow(1) & ...
    tx <= timeWindow(2);


t = ...
    tx(idx);

x = ...
    x(idx);

yd = ...
    yd(idx);


% =========================================================================
% ELIMINAR NaN / INF
% =========================================================================

valid = ...
    isfinite(t) & ...
    isfinite(x) & ...
    isfinite(yd);


t = ...
    t(valid);

x = ...
    x(valid);

yd = ...
    yd(valid);


if numel(t) < 2

    error( ...
        'No existen suficientes muestras para calcular fase.');

end


% =========================================================================
% ELIMINAR DC
% =========================================================================

x = ...
    x - mean(x);

yd = ...
    yd - mean(yd);


% =========================================================================
% FRECUENCIA
% =========================================================================

w0 = ...
    2*pi*f0;


% =========================================================================
% PROYECCIÓN COMPLEJA
% =========================================================================

E = ...
    exp(-1j*w0*t);


X = ...
    sum(x .* E);

Y = ...
    sum(yd .* E);


if abs(X) < 1e-15

    phi_deg = NaN;

    gain_ratio = NaN;

    delay_s = NaN;

    delay_samples = NaN;

    return;

end


% =========================================================================
% FUNCIÓN COMPLEJA
% =========================================================================

H = ...
    Y/X;


% =========================================================================
% FASE
% =========================================================================

phi_deg = ...
    rad2deg(angle(H));


phi_deg = ...
    mod(phi_deg + 180,360) - 180;


% =========================================================================
% GANANCIA
% =========================================================================

gain_ratio = ...
    abs(H);


% =========================================================================
% RETARDO
% =========================================================================

phi_rad = ...
    deg2rad(phi_deg);


delay_s = ...
    -phi_rad/w0;


% =========================================================================
% RETARDO EN MUESTRAS
% =========================================================================

Ts = ...
    mean(diff(t));


delay_samples = ...
    delay_s/Ts;


end