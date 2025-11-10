%%configuración
clear;
clc;
f_carrier = 20e3;
frecuency_simulation = 10*f_carrier;
sample_time = 1/(frecuency_simulation);

%%%Parámetros físicos del motor
Step_angle = 1.8; %%pasos en grados del motor
N_phases = 2; %%numero de phases
N_steps = 360/Step_angle; %%numero de pasos del motor
N_teeths = N_steps/N_phases; %%Cantidad de dientes del rotor
P = N_teeths/2;%%Número de pares de polos

I_nom = 1; %%Corriente nominal del torque 
Tdm = 0;%%18/1000;%%Torque para que no se mueva el rotor
Thold =  330/1000; %%Torque máximo para mantener la posición
Psi = Thold/(P*I_nom);
Kt = P*Psi;


R =  3.4; %% [Ohms]
L = 6/1000; %% [L]
B = 5/1000; %%roce
J = 4.7/1000000; %%inercia del motor


Vdc = 24;



%%Ganancias KPI corriente
shi_corriente = 0.707;
tau_corriente = L/R;
wn_corriente = 10/tau_corriente;
Kp_corriente = 2*shi_corriente*wn_corriente*L-R;
Ki_corriente = (wn_corriente^2) * L;


%%Ganancias KPI velocidad
shi_w = 0.707;

wn_w = wn_corriente/20;
Kp_w = (2*shi_w*wn_w*J)/Kt;
Ki_w = (wn_w^2)*J/Kt;

