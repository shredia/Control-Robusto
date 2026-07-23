Time_simulation = 0.5;
t_ref = [0.0];
w_ref = [10];

Wm_ref = timeseries(w_ref, t_ref);
Wm_ref = setinterpmethod(Wm_ref,'zoh');

% =========================
% Torque de carga
% =========================
t_load = [0.0 0.15];
Tl_ref = [0.0 0.25 ];

Tl = timeseries(Tl_ref, t_load);
Tl = setinterpmethod(Tl, 'linear');
%%%Parámetros físicos del motor
Step_angle = 1.8; %%pasos en grados del motor
N_phases = 2; %%numero de phases
N_steps = 360/Step_angle; %%numero de pasos del motor
N_teeths = N_steps/N_phases; %%Cantidad de dientes del rotor
P = N_teeths/2;%%Número de pares de polos

I_nom = 1.7; %%Corriente nominal del torque 
I_max = sqrt(2)*I_nom;
Tdm = 22/1000;%%22/1000;%%Torque para que no se mueva el rotor
Thold =  392/1000; %%Torque máximo para mantener la posición
% Cálculo de Kt basado en 2 fases excitadas
Kt = Thold / (sqrt(2) * I_nom); 
Psi = Kt/P;
Ke = Kt;
max_step_rate = 3000;
R =  2.5; %% [Ohms]
L = 6/1000; %% [L]
Phi = pi/2;
    
%%Saliencia
Lq = 9.61/1000;
Ld = 3.66/1000;
L0 = 6.635/1000;
L2 = 2.975/1000;

% Ld = L;
% Lq = Ld;
% L2 = 0;
% L0 = L;
%%Propiedades del motor (variables que creo conocer)
B_internal = 1e-4; %%roce
J_internal = 54e-7; %%inercia del motor 
%%Propiedades de carga del motor (variables que desconozco)
J_external = J_internal*0;
B_external = B_internal*0;
%%Variables del modelo (Variables que asumo más grandes) 
J_var = J_internal;
B_var = B_internal;

B_real = B_internal + B_external;
J_real = J_internal + J_external;

J_rate = J_var/J_real;
B_rate = B_var/B_real;
%%Variables de roce que desconozco
Tc = 0.002; %%Coulomb friction
Tba = 0;2.5*Tc; %%Breakwat friction
Wba = 0;2;%%breakway frictions (rad/s)


InitialSpeed = 0; %%rad/s
RI_flag = 0;





Vdc = 24;

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
Ts_Wm = Ts_current*10;
Ts_pos = Ts_Wm*10;
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

shi_w = 1;
wn_w = 2*pi*fbw_Wm;

Kp_w = (2*shi_w*wn_w*J_var )/Kt;
Ki_w = (wn_w^2)*J_var/Kt;
Kd_w = Kp_w/1000;
Tf_w = 10/1000;


Kp_pos = 1.0;    % Ganancia Proporcional (rigidez)
Ki_pos = 0;     % Ganancia Integral (elimina error en estado estable)
Kd_pos = 0.1;    % Ganancia Derivativa (amortiguamiento ante cogging)


%% ==========================
% Resonant Controller
%% ==========================


zeta_p = 0.01;
zeta_z = 0.9;

z0 = 0.98;
z6 = 0.7;

Kri = 0.001;



wr = P*w_ref;

wp = wr/sqrt(1-2*zeta_p^2);
wz = wp;

a_ri = 2*exp(-Ts_Wm*zeta_z*wz)* ...
       cos(Ts_Wm*wz*sqrt(1-zeta_z^2));

b_ri = exp(-2*Ts_Wm*zeta_z*wz);

c_ri = 2*exp(-Ts_Wm*zeta_p*wp)* ...
       cos(Ts_Wm*wp*sqrt(1-zeta_p^2));

d_ri = exp(-2*Ts_Wm*zeta_p*wp);

Kr_norm = Kri*(1-c_ri+d_ri)/(1-a_ri+b_ri);


fprintf('wr      = %.6f\n',wr);
fprintf('a_ri    = %.6f\n',a_ri);
fprintf('b_ri    = %.6f\n',b_ri);
fprintf('c_ri    = %.6f\n',c_ri);
fprintf('d_ri    = %.6f\n',d_ri);

Kr_norm = Kri*(1-c_ri+d_ri)/(1-a_ri+b_ri);

fprintf('Kr_norm = %.6f\n',Kr_norm);

%% =========================================================================
% SINTONIZACIÓN INDEPENDIENTE DEL OBSERVADOR DE PERTURBACIONES (DOB)
%% =========================================================================
fbw_DO = 30;                  % Hz (Ancho de banda del DOB, subido de ~10Hz a 120Hz)
wn_do = 2*pi*fbw_DO;           % rad/s

% Usamos una sintonización de polos repetidos en -wn_do (Polinomio de Hurwitz)
% para garantizar estabilidad estricta y convergencia rápida.
BJ_hat = B_var / J_var;

B0 = 3*wn_do - BJ_hat;
B1 = 3*wn_do^2 - B0*BJ_hat;
B2 = (wn_do^3)*(J_var/P);        % Ajuste de ganancia pura de perturbación

% Si usas la estructura exacta del paper (unidades eléctricas):
% Asegúrate de pasar estos valores calculados (beta0_calc, beta1_calc, beta2_calc)
% directo a los puertos de entrada de tu S-Function.

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
params.Nr = N_steps;

%% Cálculo de parámetros del HFI
HFI_flag = 1;
Amplitud_HFI = 0.5;      % Voltios
f_h = 1000;             % Frecuencia de inyección en Hz
wh = 2 * pi * f_h;     % Frecuencia en rad/s (para la demodulación)

% Parámetros del Filtro BPF
Fs_control = f_ekf;    % Asumiendo frecuencia de muestreo de 10kHz (ajustar si es otra)
BW_hfi = 100;          % Ancho de banda del filtro en Hz (típicamente 100-200Hz)

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
% 
% 
% %% ====================================================
% % PI convencional
% %% ====================================================
% 
% s = tf('s');
% 
% % Planta mecánica
% Gp = 1/(J_var*s + B_var);
% 
% % PI velocidad
% C_PI = Kp_w + Ki_w/s;
% 
% % Lazo abierto
% L_PI = Kt*C_PI*Gp;
% 
% % Perturbación torque -> velocidad
% Gd_PI = -Gp/(1 + L_PI);
% 
% figure
% bode(Gd_PI)
% grid on
% title('Perturbación Torque -> Velocidad (PI)')
% 
% %% ====================================================
% % Resonador
% %% ====================================================
% 
% wr = P*w_ref;
% 
% zeta_p = 0.01;
% zeta_z = 0.9;
% 
% Kri = 0.03;
% 
% R1 = Kri * ...
%     (s^2 + 2*zeta_z*wr*s + wr^2) / ...
%     (s^2 + 2*zeta_p*wr*s + wr^2);
% 
% %% ====================================================
% % Controlador total
% %% ====================================================
% 
% C_RI = C_PI + R1;
% 
% L_RI = Kt*C_RI*Gp;
% 
% Gd_RI = -Gp/(1 + L_RI);
% 
% figure
% 
% bode(Gd_PI,'r',Gd_RI,'b')
% 
% grid on
% 
% legend('PI','PI + RI')
% 
% title('Rechazo perturbación torque')
% 
% wd = P*w_ref;
% 
% [mag_PI,~] = bode(Gd_PI,wd);
% [mag_RI,~] = bode(Gd_RI,wd);
% 
% mag_PI = squeeze(mag_PI);
% mag_RI = squeeze(mag_RI);
% 
% fprintf('\n');
% fprintf('===== COMPARACION =====\n');
% 
% fprintf('PI     : %.2f dB\n',20*log10(mag_PI));
% fprintf('PI+RI  : %.2f dB\n',20*log10(mag_RI));
% 
% fprintf('Mejora : %.2f dB\n', ...
%     20*log10(mag_RI/mag_PI));