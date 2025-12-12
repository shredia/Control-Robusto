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
Tdm = 17/1000;%%18/1000;%%Torque para que no se mueva el rotor
Thold =  330/1000; %%Torque máximo para mantener la posición
Psi = Thold/(P*I_nom);
Kt = P*Psi;
Ke = Psi;

R =  3.4; %% [Ohms]
L = 6/1000; %% [L]
B = 5/10000; %%roce
J = 4.7/1000000; %%inercia del motor


Vdc = 24;

% f_carrier = 20e3;
% frecuency_simulation = 40e3;
fp = round(R/(2*pi*L)/10)*10 
fbw = 500
wd = round(fbw*2*pi/100)*100
f_carrier = 20e3
frecuency_simulation = 40e3
sample_time = 1/(frecuency_simulation);
Ts_ekf2 = 1e-3;

%%Ganancias PID Id


Kp_d = wd*L;
Ki_d = wd*R;


%%Ganancias KPI corriente


Kp_q = wd*L;
Ki_q = wd*R;


%%Ganancias KPI velocidad

shi_w = 0.3;
K_w = Kt/B;
tau_w = J/B;
wn_w = wd/10
Kp_w = (2*shi_w*wn_w*tau_w*J)/Kt;
Ki_w = (wn_w^2)*J/Kt;


B0 = 3*wn_w -(B/J);
B1 = 3*wn_w^2 -(B0*B/J);
B2 = -(J/P)*(wn_w^3);