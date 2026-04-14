clear all; clc;

%% 1. Parámetros del Motor (Tus datos)
R = 3.4;
Ld = 8.71e-3;   % Henrios
Lq = 3.18e-3;   % Henrios
P = 50;         % Pares de polos
I_nom = 1; %%Corriente nominal del torque 
Tdm = 18/1000;%%18/1000;%%Torque para que no se mueva el rotor
Thold =  330/1000; %%Torque máximo para mantener la posición
Psi = Thold/(P*I_nom);
Kt = P*Psi;
Ke = Psi;    % V/(rad/s)
J = 4.7e-6;     % Inercia
B = 5e-4;       % Fricción

% Inductancias medias y de saliencia
L0 = (Ld + Lq)/2;
L2 = (Ld - Lq)/2;

%% 2. Configuración de la Simulación
tspan = [0 0.1];             % 100ms de simulación
x0 = [0; 0; 0; 0];           % [Ia; Ib; Wm; Theta_m]
V_amp = 12;                  % Voltaje aplicado
f_Hz = 10;                   % Frecuencia de la fuente (Hz)

%% 3. Resolución de la EDO
% Definimos la función del motor
motor_ode = @(t, x) pmsm_2ph_dynamics(t, x, R, L0, L2, P, Ke, J, B, V_amp, f_Hz);

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
[t, x] = ode45(motor_ode, tspan, x0, options);

%% 4. Gráficos
figure;
subplot(3,1,1); plot(t, x(:,1), t, x(:,2)); title('Corrientes Ia, Ib'); ylabel('Amperes');
subplot(3,1,2); plot(t, x(:,3)); title('Velocidad Mecánica \omega_m'); ylabel('rad/s');
subplot(3,1,3); plot(t, mod(x(:,4)*P, 2*pi)); title('Ángulo Eléctrico \theta_e'); ylabel('rad');

%% --- Función de Dinámica ---
function dxdt = pmsm_2ph_dynamics(t, x, R, L0, L2, P, Ke, J, B, V_amp, f_Hz)
    % Estados
    Ia = x(1); Ib = x(2); Wm = x(3); Th_m = x(4);
    Th_e = P * Th_m;
    We = P * Wm;

    % Voltajes de entrada (Simulamos una fuente senoidal simple)
    Va = V_amp * cos(2*pi*f_Hz*t);
    Vb = V_amp * sin(2*pi*f_Hz*t);

    % 1. Matriz de Inductancia Variable L(theta)
    La  = L0 + L2*cos(2*Th_e);
    Lb  = L0 - L2*cos(2*Th_e);
    Lab = L2*sin(2*Th_e);
    
    detL = La*Lb - Lab^2; % Debe ser Ld*Lq (siempre positivo)

    % 2. Derivadas de Inductancia (dL/dt = dL/dth * dth/dt)
    dLa_dt  = -2*P*Wm * L2*sin(2*Th_e);
    dLb_dt  =  2*P*Wm * L2*sin(2*Th_e);
    dLab_dt =  2*P*Wm * L2*cos(2*Th_e);

    % 3. Back-EMF (Fuerza Contraelectromotriz)
    Ea = -Ke * Wm * sin(Th_e);
    Eb =  Ke * Wm * cos(Th_e);

    % 4. Voltajes Netos (Sin incluir dI/dt)
    % V_net = V - R*I - E - I*(dL/dt)
    Va_net = Va - R*Ia - Ea - (Ia*dLa_dt + Ib*dLab_dt);
    Vb_net = Vb - R*Ib - Eb - (Ib*dLb_dt + Ia*dLab_dt);

    % 5. Ecuaciones Eléctricas (Resolviendo el sistema acoplado)
    dIa = ( Lb*Va_net - Lab*Vb_net) / detL;
    dIb = ( La*Vb_net - Lab*Va_net) / detL;

    % 6. Torque Eléctrico (Imán + Reluctancia)
    Te = P * ( (Ke/P)*(Ib*cos(Th_e) - Ia*sin(Th_e)) + ...
         (L0-L2 - (L0+L2))*( (Ib^2 - Ia^2)*sin(2*Th_e)/2 + Ia*Ib*cos(2*Th_e) ) );
     % Nota: Simplifiqué (Ld-Lq) usando tus variables L0 y L2

    % 7. Ecuaciones Mecánicas
    dWm = (Te - B*Wm) / J;
    dTh = Wm;

    dxdt = [dIa; dIb; dWm; dTh];
end