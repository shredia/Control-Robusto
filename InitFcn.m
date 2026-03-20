Time_simulation = 0.5;
t_ref = [0.0];
w_ref = [45];

Wm_ref = timeseries(w_ref, t_ref);
Wm_ref = setinterpmethod(Wm_ref,'zoh');
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
%%Propiedades del motor (variables que creo conocer)
B_internal = 5/10000; %%roce
J_internal = 4.7/1000000; %%inercia del motor 
%%Propiedades de carga del motor (variables que desconozco)
J_external = J_internal*1;
B_external = B_internal*1;
%%Variables del modelo (Variables que asumo más grandes) 
J_var = J_internal*2;
B_var = B_internal*2;

B_real = B_internal + B_external;
J_real = J_internal + J_external;

J_rate = J_var/J_real;
B_rate = B_var/B_real;
%%Variables de roce que desconozco
Tc = 0.002; %%Coulomb friction
Tba = 2.5*Tc; %%Breakwat friction
Wba = 2;%%breakway frictions (rad/s)

InitialSpeed = 0; %%rad/s




Vdc = 24;

fbw = 800;          % Hz (BW corriente)
wd  = 2*pi*fbw;      % rad/s  <-- CORRECTO

f_carrier = 40e3;
frecuency_simulation = 100e3;
sample_time = 1/frecuency_simulation;


f_current = f_carrier;
f_Wm = 30;
f_ekf = f_current;

Ts_ekf = 1/f_ekf;


Ts_current = 1/f_current;
Ts_Wm = Ts_current*10;
Ts_DO = Ts_Wm;
Kp_d = wd*L;
Ki_d = wd*R;

Kp_q = wd*L;
Ki_q = wd*R;

shi_w = 4;
wn_w = 2*pi*f_Wm

Kp_w = (2*shi_w*wn_w*J_var)/Kt;
Ki_w = (wn_w^2)*J_var/Kt;

% Si B desconocido, no uses B0,B1 con (B/J) real
BJ_hat = B_var/J_var;     % o una estimación
B0 = 3*wn_w - BJ_hat;
B1 = 3*wn_w^2 - B0*BJ_hat;
B2 = -(J_var/P)*(wn_w^3);