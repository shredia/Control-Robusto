%% ============================================================
% LUT ANTICOGGING:
% theta_m mecánico  -->  Iq_cogging
%% ============================================================


% ---------------- CONFIGURACIÓN ----------------
N_LUT       = 720;   % puntos por vuelta mecánica
t_ini       = 2.0;   % [s] inicio de régimen permanente
theta_offset = 0;    % offset angular, inicialmente cero
smooth_bins = 9;     % debe ser impar
signo_LUT   = +1;    % cambiar a -1 si la compensación empeora el rizado

% ---------------- EXTRAER DATOS ----------------
theta = squeeze(theta_M_planta_M2.Data);
t_theta = theta_M_planta_M2.Time(:);

iq_cog = squeeze(Iq_M2.Data);
t_iq = Iq_M2.Time(:);

theta = theta(:);
iq_cog = iq_cog(:);

% Llevar Iq_cogging al mismo vector temporal de theta_m
if length(t_iq) ~= length(t_theta) || any(t_iq ~= t_theta)
    iq_cog = interp1(t_iq, iq_cog, t_theta, 'linear', 'extrap');
end

% ---------------- ELIMINAR TRANSITORIO ----------------
idx_ss = t_theta >= t_ini & isfinite(theta) & isfinite(iq_cog);

t_ss     = t_theta(idx_ss);
theta_ss = theta(idx_ss);
iq_ss    = signo_LUT * iq_cog(idx_ss);

% Posición dentro de una revolución mecánica
theta_mod = mod(theta_ss - theta_offset, 2*pi);

% ---------------- PROMEDIO POR SECTOR ANGULAR ----------------
Iq_LUT = nan(N_LUT,1);
N_muestras = zeros(N_LUT,1);

% Cada posición k representa una zona de ángulo mecánico
idx_bin = floor(theta_mod/(2*pi)*N_LUT) + 1;
idx_bin(idx_bin > N_LUT) = N_LUT;

for k = 1:N_LUT

    datos_k = iq_ss(idx_bin == k);

    if ~isempty(datos_k)
        % La mediana reduce el efecto de picos o ruido
        Iq_LUT(k) = median(datos_k);
        N_muestras(k) = length(datos_k);
    end
end

% ---------------- RELLENAR ZONAS SIN DATOS ----------------
idx_valid = ~isnan(Iq_LUT);

if sum(idx_valid) < 2
    error('No hay suficientes datos para crear la LUT. Aumenta Time_simulation.');
end

x = (1:N_LUT).';
x_valid = x(idx_valid);
y_valid = Iq_LUT(idx_valid);

% Interpolación circular: conecta correctamente 0 con 2*pi
x_ext = [x_valid - N_LUT; x_valid; x_valid + N_LUT];
y_ext = [y_valid; y_valid; y_valid];

Iq_LUT = interp1(x_ext, y_ext, x, 'linear');

% ---------------- SUAVIZADO CIRCULAR ----------------
if mod(smooth_bins,2) == 0
    error('smooth_bins debe ser impar.');
end

n_extra = floor(smooth_bins/2);

Iq_pad = [Iq_LUT(end-n_extra+1:end);
          Iq_LUT;
          Iq_LUT(1:n_extra)];

Iq_pad = movmean(Iq_pad, smooth_bins);

Iq_LUT = Iq_pad(n_extra+1:n_extra+N_LUT);

% Eliminar componente DC: la LUT no debe agregar torque promedio
Iq_LUT = Iq_LUT - mean(Iq_LUT);

% ---------------- TABLA PARA SIMULINK ----------------
Breakpoints_theta = linspace(0, 2*pi, N_LUT+1).';

% Se repite el primer valor al final para cerrar la periodicidad
TableData_Iq_comp = [Iq_LUT; Iq_LUT(1)];

save('LUT_Iq_cogging.mat', ...
    'Breakpoints_theta', ...
    'TableData_Iq_comp', ...
    'theta_offset', ...
    'N_muestras');

% ---------------- GRÁFICOS ----------------
figure

subplot(2,1,1)
plot(Breakpoints_theta, TableData_Iq_comp, 'LineWidth', 1.2)
grid on
xlabel('\theta_m mod 2\pi [rad]')
ylabel('I_{q,comp} [A]')
title('LUT de compensación de cogging')

subplot(2,1,2)
plot(Breakpoints_theta(1:end-1), N_muestras, 'LineWidth', 1.2)
grid on
xlabel('\theta_m mod 2\pi [rad]')
ylabel('Muestras por sector')
title('Cobertura angular de datos')