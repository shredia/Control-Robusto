%% ================================================================
% ANALISIS DE FASE Wm REAL vs EKF
% ================================================================

Nc = 200;

t_ini = 0.04;
t_fin = 0.10;

%% ---------------------------------------------------------------
% Señales originales
% NO eliminar la media aquí
% ---------------------------------------------------------------

Wm_real = Wm_planta_M2;
Wm_obs  = Wm_obs_M2;

%% ---------------------------------------------------------------
% Extraer datos de timeseries
% ---------------------------------------------------------------

tr = Wm_real.Time(:);
Wr = squeeze(Wm_real.Data);
Wr = Wr(:);

to = Wm_obs.Time(:);
Wo = squeeze(Wm_obs.Data);
Wo = Wo(:);

%% ---------------------------------------------------------------
% Seleccionar ventana
% ---------------------------------------------------------------

idx_r = (tr >= t_ini) & (tr <= t_fin);
idx_o = (to >= t_ini) & (to <= t_fin);

tr = tr(idx_r);
Wr = Wr(idx_r);

to = to(idx_o);
Wo = Wo(idx_o);

%% ---------------------------------------------------------------
% Interpolar si los vectores temporales no coinciden
% ---------------------------------------------------------------

if length(tr) ~= length(to) || any(abs(tr-to) > 1e-12)

    Wo = interp1(to,Wo,tr,'linear');

end

t = tr;

%% ---------------------------------------------------------------
% Velocidad media REAL
% ---------------------------------------------------------------

Wr_mean = mean(Wr);
Wo_mean = mean(Wo);

fprintf('\nWm real media = %.4f rad/s\n',Wr_mean);
fprintf('Wm EKF media  = %.4f rad/s\n',Wo_mean);

%% ---------------------------------------------------------------
% Ripple de velocidad
%
% SOLO AHORA quitamos la componente DC
% ---------------------------------------------------------------

wr = Wr - Wr_mean;
wo = Wo - Wo_mean;

%% ---------------------------------------------------------------
% Correlacion sin toolbox
% ---------------------------------------------------------------

rho = sum(wr .* wo) / ...
      sqrt(sum(wr.^2) * sum(wo.^2));

fprintf('Correlacion ripple = %.4f\n',rho);

%% ---------------------------------------------------------------
% Frecuencia esperada de cogging
% ---------------------------------------------------------------

wcog = Nc * abs(Wr_mean);

fcog = wcog/(2*pi);

fprintf('\nwcog teorica = %.2f rad/s\n',wcog);
fprintf('fcog teorica = %.2f Hz\n',fcog);

%% ---------------------------------------------------------------
% Fase directamente a la frecuencia de cogging
% ---------------------------------------------------------------

E = exp(-1j*2*pi*fcog*t);

Xreal = sum(wr .* E);
Xobs  = sum(wo .* E);

phi_real = angle(Xreal);
phi_obs  = angle(Xobs);

dphi = angle(exp(1j*(phi_obs - phi_real)));

fprintf('\nFase real = %.2f grados\n',phi_real*180/pi);
fprintf('Fase EKF  = %.2f grados\n',phi_obs*180/pi);

fprintf('\n========================================\n');
fprintf('DESFASE EKF - REAL = %.2f grados\n', ...
        dphi*180/pi);
fprintf('========================================\n');


%% ================================================================
% BUSCAR RETARDO POSITIVO DEL EKF
% ================================================================

Ts = mean(diff(t));

maxLag = 25;

lags = 0:maxLag;
Rlag = zeros(size(lags));

for k = 1:length(lags)

    L = lags(k);

    if L == 0

        a = wr;
        b = wo;

    else

        % Wobs ocurre L muestras después que Wreal
        a = wr(1:end-L);
        b = wo(1+L:end);

    end

    Rlag(k) = sum(a.*b) / ...
        sqrt(sum(a.^2)*sum(b.^2));

end

[Rmax,idx] = max(Rlag);

lag_opt = lags(idx);

delay_opt = lag_opt*Ts;

phi_opt = -360*fcog*delay_opt;

fprintf('\n========================================\n');
fprintf('Mejor correlacion causal = %.4f\n',Rmax);
fprintf('Retardo EKF = %d samples\n',lag_opt);
fprintf('Retardo EKF = %.4f ms\n',delay_opt*1000);
fprintf('Fase equivalente = %.2f grados\n',phi_opt);
fprintf('========================================\n');