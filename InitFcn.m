%%configuración
clearvars;
clear;
clc;

%%%Parámetros físicos del motor
Step_angle = 1.8; %%pasos en grados del motor
N_phases = 2; %%numero de phases
N_steps = 360/Step_angle; %%numero de pasos del motor
N_teeths = N_steps/N_phases; %%Cantidad de dientes del rotor
P = N_teeths/2;%%Número de pares de polos

I_nom = 1; %%Corriente nominal del torque 
Tdm = 18/1000;%%18/1000;%%Torque para que no se mueva el rotor
Thold =  330/1000; %%Torque máximo para mantener la posición
Psi = Thold/(P*I_nom);
Kt = P*Psi;
Ke = Psi;
max_step_rate = 3000;
R =  3.4; %% [Ohms]
L = 6/1000; %% [L]
B_internal = 5/10000; %%roce
J_internal = 4.7/1000000; %%inercia del motor 
J_external = J_internal*0.1;
B_external = B_internal*0.2;
Tc = 0.002; %%Coulomb friction
Tba = 2.5*Tc; %%Breakwat friction
Wba = 2;%%breakway frictions (rad/s)

InitialSpeed = 0; %%rad/s




Vdc = 24;

fbw = 500;          % Hz (BW corriente)
wd  = 2*pi*fbw;      % rad/s  <-- CORRECTO

f_carrier = 40e3;
frecuency_simulation = 100e3;
sample_time = 1/frecuency_simulation;


f_current = fbw*10;
f_Wm = 30;
f_ekf = f_current;

Ts_ekf = 1/f_ekf;


Ts_current = 1/f_current;
Ts_Wm = 1/(500);
Ts_DO = Ts_Wm;
Kp_d = wd*L;
Ki_d = wd*R;

Kp_q = wd*L;
Ki_q = wd*R;

shi_w = 1.5;
wn_w = 2*pi*f_Wm

Kp_w = (2*shi_w*wn_w*J_internal)/Kt;
Ki_w = (wn_w^2)*J_internal/Kt;

% Si B desconocido, no uses B0,B1 con (B/J) real
BJ_hat = 0;     % o una estimación
B0 = 3*wn_w - BJ_hat;
B1 = 3*wn_w^2 - B0*BJ_hat;
B2 = -(J_internal/P)*(wn_w^3);