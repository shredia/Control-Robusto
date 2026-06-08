Time_simulation = 5;
t_ref = [0.0];
w_ref = [0.2];

Wm_ref = timeseries(w_ref, t_ref);
Wm_ref = setinterpmethod(Wm_ref,'zoh');
%%%Parámetros físicos del motor
Step_angle = 1.8; %%pasos en grados del motor
N_phases = 2; %%numero de phases
N_steps = 360/Step_angle; %%numero de pasos del motor
N_teeths = N_steps/N_phases; %%Cantidad de dientes del rotor
P = N_teeths/2;%%Número de pares de polos

I_nom = 1; %%Corriente nominal del torque 
I_max = sqrt(2*I_nom);
Tdm = 18/1000;%%18/1000;%%Torque para que no se mueva el rotor
Thold =  330/1000; %%Torque máximo para mantener la posición
% Cálculo de Kt basado en 2 fases excitadas
Kt = Thold / (sqrt(2) * I_nom); 
Psi = Kt/P;
Ke = Kt;
max_step_rate = 3000;
R =  2.5; %% [Ohms]
L = 6/1000; %% [L]
Phi = pi/2;
    
%%Saliencia
Ld = 8.71/1000;
Lq = 3.18/1000;
L0 = 5.945/1000;
L2 = 2.765/1000;

% Ld = L;
% Lq = Ld;
% L2 = 0;
% L0 = L;
%%Propiedades del motor (variables que creo conocer)
B_internal = 1e-4; %%roce
J_internal = 4.7/1000000; %%inercia del motor 
%%Propiedades de carga del motor (variables que desconozco)
J_external = J_internal*1;
B_external = B_internal*1;
%%Variables del modelo (Variables que asumo más grandes) 
J_var = J_internal;
B_var = B_internal;

B_real = B_internal + B_external;
J_real = J_internal + J_external;

J_rate = J_var/J_real;
B_rate = B_var/B_real;
%%Variables de roce que desconozco
Tc = 0.002; %%Coulomb friction
Tba = 2.5*Tc; %%Breakwat friction
Wba = 2;%%breakway frictions (rad/s)

InitialSpeed = 0; %%rad/s




Vdc = 12;

fbw_d = 500;          % Hz (BW corriente)
fbw_q = 500;          % Hz (BW corriente)

wd_d = 2*pi*fbw_d;      % rad/s
wd_q = 2*pi*fbw_q;      % rad/s 

fbw_Wm = 20;         % Hz (Bw Wm)
frecuency_simulation = 10e3;
f_carrier = 20e3;
sample_time = 1/frecuency_simulation;


f_current = f_carrier/2;

f_ekf = f_current;

Ts_ekf = 1/f_ekf;


Ts_current = 1/f_current;
Ts_Wm = Ts_current;
f_wm = 1/Ts_Wm;
Ts_DO = Ts_Wm;
Kp_d = wd_d*L;
Ki_d = wd_d*R;

Kp_q = wd_q*L;
Ki_q = wd_q*R;

shi_d = 0.707;
shi_q = 0.707;

Kp_q_salient = 2*shi_q*wd_q*Lq-R
Ki_q_salient = wd_q^2*Lq

Kp_d_salient = 2*shi_d*wd_d*Ld -R
Ki_d_salient = wd_d^2*Ld

shi_w = 0.707;
wn_w = 2*pi*fbw_Wm;

Kp_w = (2*shi_w*wn_w*J_var - B_var)/Kt;
Ki_w = (wn_w^2)*J_var/Kt;

% Si B desconocido, no uses B0,B1 con (B/J) real
BJ_hat = B_var/J_var;     % o una estimación
B0 = 3*wn_w - BJ_hat;
B1 = 3*wn_w^2 - B0*BJ_hat;
B2 = -(J_var/P)*(wn_w^3);

% --- Empaquetado para S-Function Level-2 ---
params.R      = R;
params.Ld     = Ld;
params.Lq     = Lq;
params.P      = P;
params.Ke     = Ke;
params.J_real = J_real; % Usamos el valor real para la planta
params.B_real = B_real; % Usamos el valor real para la planta
params.Kt     = Kt;
params.Tdm = Tdm;
params.N_steps = N_steps;
params.Phi = Phi;


%% Cálculo de parámetros del HFI
Amplitud_HFI = 0.5;      % Voltios
f_h = 1000;             % Frecuencia de inyección en Hz
wh = 2 * pi * f_h;     % Frecuencia en rad/s (para la demodulación)

% Parámetros del Filtro BPF
Fs_control = f_ekf;    % Asumiendo frecuencia de muestreo de 10kHz (ajustar si es otra)
BW_hfi = 200;          % Ancho de banda del filtro en Hz (típicamente 100-200Hz)

% Diseño del filtro de 1er orden de butterworth (produce un BPF de 2do orden)
% Normalización de frecuencias [f_inferior, f_superior] respecto a Fs/2
Wn = [(f_h - BW_hfi/2)/(Fs_control/2), (f_h + BW_hfi/2)/(Fs_control/2)];

[b_hfi, a_hfi] = butter(1, Wn, 'bandpass');

% Coeficientes para tu función EKF
b0_hfi = b_hfi(1);
b1_hfi = b_hfi(2);
b2_hfi = b_hfi(3);
a1_hfi = a_hfi(2);
a2_hfi = a_hfi(3);

fprintf('Coeficientes calculados para f_h = %d Hz:\n', f_h);
fprintf('b0 = %.6f, b1 = %.6f, b2 = %.6f\n', b0_hfi, b1_hfi, b2_hfi);
fprintf('a1 = %.6f, a2 = %.6f\n', a1_hfi, a2_hfi);

% Análisis de lazo cerrado de velocidad en Matlab
% s = tf('s');
% G_planta = Kt / (J_real * s + B_real);
% C_pi = Kp_w + Ki_w/s;
% L_a = C_pi * G_planta; % Lazo abierto
% T = feedback(L_a, 1);  % Lazo cerrado
% 
% figure
% margin(L_a)
% grid on
% title('Bode con márgenes')
% 
% figure
% step(T)
% grid on
% title('Respuesta al escalón')
% 
% figure
% pzmap(T)
% grid on
% title('Mapa polo-cero')


% %% 1. Variable de Laplace
% 
% s = tf('s');
% 
% %% 2. Frecuencia de la perturbación
% 
% wd = 4*P*w_ref;      % rad/s
% fd = wd/(2*pi);      % Hz
% 
% fprintf('Frecuencia perturbación: %.4f rad/s = %.4f Hz\n', wd, fd);
% 
% %% 3. Controlador PI de velocidad
% 
% Cw = Kp_w + Ki_w/s;
% 
% %% 4. Transferencia perturbación -> velocidad sin PI
% 
% Gd_w_open = -1/(J_var*s + B_var);
% 
% %% 5. Transferencia perturbación -> velocidad con PI
% 
% Gd_w_closed = -1/(J_var*s + B_var + Kt*Cw);
% 
% %% 6. Bode comparativo
% 
% figure;
% bode(Gd_w_open, Gd_w_closed);
% grid on;
% legend('Sin PI', 'Con PI');
% title('Bode: perturbación de torque T_d \rightarrow velocidad \omega_m');
% 
% %% 7. Marcar frecuencia de perturbación en el Bode
% 
% figure;
% bode(Gd_w_closed);
% grid on;
% title('Bode: T_d \rightarrow \omega_m con PI de velocidad');
% 
% hold on;
% 
% % Obtener magnitud y fase exactamente en wd
% [mag, phase] = bode(Gd_w_closed, wd);
% 
% mag = squeeze(mag);
% phase = squeeze(phase);
% 
% mag_db = 20*log10(mag);
% 
% fprintf('\n===== Evaluación en frecuencia de perturbación =====\n');
% fprintf('|Gd_w(jwd)| = %.6e rad/s por N*m\n', mag);
% fprintf('|Gd_w(jwd)| = %.2f dB\n', mag_db);
% fprintf('Fase = %.2f grados\n', phase);
% 
% %% 8. Estimar rizado desde el Bode
% 
% A_w = mag*Tdm;             % amplitud del rizado [rad/s]
% ripple_pp = 2*A_w;         % peak-to-peak [rad/s]
% ripple_percent = ripple_pp/w_ref*100;
% 
% fprintf('\n===== Rizado estimado =====\n');
% fprintf('A_w = %.6f rad/s\n', A_w);
% fprintf('Ripple pp = %.6f rad/s\n', ripple_pp);
% fprintf('Ripple percent = %.2f %%\n', ripple_percent);