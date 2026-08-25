% =========================================================================
% PLANTA MECÁNICA
% =========================================================================

J = motor_params.J_internal;
B = motor_params.B_internal;

Gm = tf(1, [J B]);

% Frecuencia de cogging
f_cog = 312.735;             % Hz
w_cog = 2*pi*f_cog;          % rad/s


% =========================================================================
% BODE
% =========================================================================

figure;

h = bodeplot(Gm);

grid on;

title('Planta mecánica W_m / T_e');


% =========================================================================
% GANANCIA Y FASE EN LA FRECUENCIA DE COGGING
% =========================================================================

[mag, phase] = bode(Gm, w_cog);

mag = squeeze(mag);
phase = squeeze(phase);

mag_dB = 20*log10(mag);

fprintf('\n');
fprintf('Frecuencia de cogging:\n');
fprintf('f_cog  = %.3f Hz\n', f_cog);
fprintf('w_cog  = %.3f rad/s\n', w_cog);

fprintf('\nPlanta mecánica en f_cog:\n');
fprintf('Magnitud = %.6f\n', mag);
fprintf('Magnitud = %.3f dB\n', mag_dB);
fprintf('Fase     = %.3f grados\n', phase);