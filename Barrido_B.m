%% Barrido de J - 3 motores (M0, M1, M2)
clearvars;
clc;

InitFcn();   % Carga parámetros base del modelo

%% =========================
%% Configuración general
%% =========================
modelName = "Simulacion_stepper";

B_factors = [B_var];   % multiplicadores respecto a B_internal
N = numel(B_factors);

Fmax      = 2000;
Kmax      = 5;
doPlotFFT = true;   % false durante barrido, para no abrir muchas figuras

t_inicial = 0.2;
t_final   = 0.4;
Time_simulation = 0.5;

motorNames = ["M2"];

%% =========================
%% Carpeta de guardado
%% =========================
rootFolder = "barrido_B";
dataFolder = fullfile(rootFolder, "data_cases");
resFolder  = fullfile(rootFolder, "res_cases");


if ~exist(rootFolder, 'dir')
    mkdir(rootFolder);
end
if ~exist(dataFolder, 'dir')
    mkdir(dataFolder);
end
if ~exist(resFolder, 'dir')
    mkdir(resFolder);
end


%% =========================
%% Referencia de velocidad
%% =========================
t_ref = [0.0];
w_ref = [1];

Wm_ref = timeseries(w_ref, t_ref);
Wm_ref = setinterpmethod(Wm_ref,'zoh');

% como Wm_ref es constante, fe_ref es simplemente:
fe_ref = mean(Wm_ref.Data) * P / (2*pi);

%% =========================
%% Prealocación
%% =========================
data = struct();
res  = struct();

summaryStruct = struct();
row = 0;

%% =========================
%% 1) LOOP DE SIMULACIONES
%% =========================
for k = 1:N

    fprintf('\n=============================\n');
    fprintf('Simulando caso %d de %d | Bx%.3g\n', k, N, B_factors(k));
    fprintf('=============================\n');

    simIn = Simulink.SimulationInput(modelName);

    J_real = J_internal + J_external;
    B_real = B_internal + B_external;

    J_var = J_real;
    B_var = B_internal* B_factors(k);   % fijo, según tu archivo original

    J_rate = J_var / J_real;
    B_rate = B_var / B_real;

    Kp_w = (2*shi_w*wn_w*J_var)/Kt;
    Ki_w = (wn_w^2)*J_var/Kt;

    % Si B desconocido, no uses B0,B1 con (B/J) real
    BJ_hat = 0.5;
    B0 = 3*wn_w - BJ_hat;
    B1 = 3*wn_w^2 - B0*BJ_hat;
    B2 = -(J_var/P)*(wn_w^3);

    % Mantengo tu enfoque con base workspace
    assignin('base','J_var',J_var);
    assignin('base','B_var',B_var);
    assignin('base','J_rate',J_rate);
    assignin('base','B_rate',B_rate);
    assignin('base','Kp_w',Kp_w);
    assignin('base','Ki_w',Ki_w);
    assignin('base','B0',B0);
    assignin('base','B1',B1);
    assignin('base','B2',B2);
    assignin('base','Wm_ref',Wm_ref);


    out = sim(simIn);

    %% ---- Guardar metadatos en memoria ----
    data(k).name      = "Bx" + string(B_factors(k));
    data(k).B_factor  = B_factors(k);
    data(k).J_var     = J_var;
    data(k).B_var     = B_var;
    data(k).J_rate    = J_rate;
    data(k).B_rate    = B_rate;
    data(k).Kp_w      = Kp_w;
    data(k).Ki_w      = Ki_w;
    data(k).B0        = B0;
    data(k).B1        = B1;
    data(k).B2        = B2;
    data(k).simOut    = out;
    data(k).time      = out.tout;

    %% ---- Guardar señales de los 3 motores ----
    % data(k).M0.Ia = out.Ia_planta_M0;
    % data(k).M0.Ib = out.Ib_planta_M0;
    % data(k).M0.We = out.We_planta_M0;
    % data(k).M0.Wm = out.Wm_planta_M0;
    % data(k).M0.Te = out.Te_planta_M0;
    % 
    % data(k).M1.Ia = out.Ia_planta_M1;
    % data(k).M1.Ib = out.Ib_planta_M1;
    % data(k).M1.We = out.We_planta_M1;
    % data(k).M1.Wm = out.Wm_planta_M1;
    % data(k).M1.Te = out.Te_planta_M1;

    data(k).M2.Ia = out.Ia_planta_M2;
    data(k).M2.Ib = out.Ib_planta_M2;
    data(k).M2.We = out.We_planta_M2;
    data(k).M2.Wm = out.Wm_planta_M2;
    data(k).M2.Te = out.Te_planta_M2;

      % ===== Señales del observador M2 =====
    % Asegúrate de que estas señales existan en tu modelo:
    data(k).M2.theta_e     = out.theta_e_planta_M2;   % ángulo real
    data(k).M2.theta_e_obs = out.theta_e_obs_M2;      % ángulo observado
    data(k).M2.Wm_obs      = out.Wm_obs_M2;           % velocidad observada

    data(k).M2.err_theta_e_deg = out.err_theta_deg;
    data(k).M2.err_Wm = out.err_Wm;

    %% ---- Guardar archivo data por simulación ----
    data_k = data(k); %#ok<NASGU>
    save(fullfile(dataFolder, "data_" + data(k).name + ".mat"), 'data_k', '-v7.3');

    %% =========================
    %% 2) ANÁLISIS FFT DEL CASO k
    %% =========================
    name_k = data(k).name;

    res(k).name     = name_k;
    res(k).B_factor = data(k).B_factor;
    res(k).B_var    = data(k).B_var;

    for m = 1:numel(motorNames)

        mot = motorNames(m);

        res(k).(mot).Ia = fft_analyzer(data(k).(mot).Ia, ...
            'Name', "Ia_" + mot + "_" + name_k, ...
            'MaxPeaks', 20, ...
            'MinPeakFrac', 0.05, ...
            'Fref', fe_ref, ...
            'Fmax', Fmax, ...
            'Kmax', Kmax, ...
            'TimeWindow', [t_inicial t_final], ...
            'DoPlot', doPlotFFT);

        res(k).(mot).Ib = fft_analyzer(data(k).(mot).Ib, ...
            'Name', "Ib_" + mot + "_" + name_k, ...
            'MaxPeaks', 20, ...
            'MinPeakFrac', 0.05, ...
            'Fref', fe_ref, ...
            'Fmax', Fmax, ...
            'Kmax', Kmax, ...
            'TimeWindow', [t_inicial t_final], ...
            'DoPlot', doPlotFFT);

        res(k).(mot).We = fft_analyzer(data(k).(mot).We, ...
            'Name', "We_" + mot + "_" + name_k, ...
            'MaxPeaks', 20, ...
            'MinPeakFrac', 0.05, ...
            'Fref', fe_ref, ...
            'Fmax', Fmax, ...
            'Kmax', Kmax, ...
            'TimeWindow', [t_inicial t_final], ...
            'DoPlot', doPlotFFT);

        res(k).(mot).Wm = fft_analyzer(data(k).(mot).Wm, ...
            'Name', "Wm_" + mot + "_" + name_k, ...
            'MaxPeaks', 20, ...
            'MinPeakFrac', 0.05, ...
            'Fref', fe_ref, ...
            'Fmax', Fmax, ...
            'Kmax', Kmax, ...
            'TimeWindow', [t_inicial t_final], ...
            'DoPlot', doPlotFFT);

        res(k).(mot).Te = fft_analyzer(data(k).(mot).Te, ...
            'Name', "Te_" + mot + "_" + name_k, ...
            'MaxPeaks', 20, ...
            'MinPeakFrac', 0.05, ...
            'Fref', fe_ref, ...
            'Fmax', Fmax, ...
            'Kmax', Kmax, ...
            'TimeWindow', [t_inicial t_final], ...
            'DoPlot', doPlotFFT);

        % frecuencia eléctrica medida desde We
        res(k).(mot).fe = res(k).(mot).We.mean / (2*pi);


        %% ---- errores de observación (solo M2) ----
error_theta_deg_mean = NaN;
error_theta_deg_rms  = NaN;
error_Wm_mean        = NaN;
error_Wm_rms         = NaN;

if mot == "M2"

    % tomar errores desde data
    error_theta_deg_mean = data(k).M2.err_theta_e_deg;
    error_Wm_mean        =  data(k).M2.err_Wm;


    % Guardar también en res
    res(k).M2.err_theta_deg     = error_theta_deg_mean;
    res(k).M2.err_Wm            = error_Wm_mean;

end

%% ---- llenar summaryStruct ----
row = row + 1;

summaryStruct(row).Modelo              = string(name_k);
summaryStruct(row).Motor               = string(mot);
summaryStruct(row).J_factor            = data(k).B_factor;
summaryStruct(row).J_var               = data(k).B_var;
summaryStruct(row).fe_Hz               = res(k).(mot).fe;
summaryStruct(row).THD_Ia_pct          = 100 * res(k).(mot).Ia.THD;
summaryStruct(row).THD_Ib_pct          = 100 * res(k).(mot).Ib.THD;
summaryStruct(row).Te_mean             = res(k).(mot).Te.mean;
summaryStruct(row).Te_ripple_rms       = res(k).(mot).Te.ripple_rms_ac;
summaryStruct(row).Te_ripple_pp        = res(k).(mot).Te.ripple_pp;
summaryStruct(row).Wm_mean             = res(k).(mot).Wm.mean;
summaryStruct(row).Wm_ripple_pp        = res(k).(mot).Wm.ripple_pp;
summaryStruct(row).error_theta_deg     = error_theta_deg_mean;
summaryStruct(row).error_Wm            = error_Wm_mean;

    end

    %% ---- Guardar archivo res por simulación ----
    res_k = res(k); %#ok<NASGU>
    save(fullfile(resFolder, "res_" + data(k).name + ".mat"), 'res_k', '-v7.3');

end

%% =========================
%% 3) TABLA RESUMEN FINAL
%% =========================
summaryTable = struct2table(summaryStruct);
disp(summaryTable)

%% =========================
%% 4) GUARDAR RESUMEN GLOBAL
%% =========================
save(fullfile(rootFolder, 'summaryStruct.mat'), 'summaryStruct');
save(fullfile(rootFolder, 'summaryTable.mat'), ...
    'summaryTable', 'B_factors', 'Wm_ref', 'fe_ref', 't_inicial', 't_final');

fprintf('\nGuardado completado en carpeta: %s\n', rootFolder);


%% =========================
%% FUNCIÓN AUXILIAR
%% =========================
function x = getTsData(sig)
    if isa(sig, 'timeseries')
        x = sig.Data;
    else
        x = sig;
    end
    x = squeeze(x);
end