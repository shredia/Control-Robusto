%% ============================================================
% ANALISIS COMPLETO DE GANANCIAS Y DESFASES
% PMSM / STEPPER SENSORLESS
%
% Evalua cada bloque en la frecuencia de cogging:
%
% PI_Wm -> MTPA -> lazo Id/Iq -> torque -> mecanica -> filtro Wm
%
% IMPORTANTE:
% El MTPA se considera algebraico. Si se desea mayor precision,
% se puede introducir posteriormente su ganancia incremental exacta.
% ============================================================


%% ============================================================
% 1. FRECUENCIA DE ANALISIS
% ============================================================

f_cog = 312.735;              % [Hz]
w_cog = 2*pi*f_cog;           % [rad/s]

s = 1j*w_cog;


%% ============================================================
% 2. PARAMETROS MECANICOS
% ============================================================

J = motor_params.J_internal;                   % [kg*m^2]
B = motor_params.B_internal;                     % [N*m*s/rad]


%% ============================================================
% 3. PARAMETROS ELECTRICOS
% ============================================================

R  = motor_params.R;                       % [ohm] <-- CAMBIAR POR TU VALOR REAL

Ld = motor_params.Ld;                 % [H]
Lq = motor_params.Lq;                 % [H]

P   = motor_params.P;                     % pares de polos
Psi = motor_params.Ke;                 % flujo / constante equivalente
                              % AJUSTAR si tu modelo usa otro valor


%% ============================================================
% 4. PI DE VELOCIDAD
% ============================================================

% ---------- PON AQUI TUS GANANCIAS REALES ----------
Kp_w = PI_MTPA_params.Kp_w;                  % <-- cambiar
Ki_w = PI_MTPA_params.Ki_w;                   % <-- cambiar
Kd_w = PI_MTPA_params.Kd_w;

% PI(s) = Kp + Ki/s
C_w = Kp_w + Ki_w/s + Kd_w*s;


%% ============================================================
% 5. FILTRO DE VELOCIDAD
% ============================================================

fbw_Wm_filter = PI_MTPA_params.fbw_Wm;           % [Hz]
w_filter = 2*pi*fbw_Wm_filter;

% Filtro primer orden
%
% H(s) = wf / (s + wf)

H_Wm = w_filter/(s + w_filter);


%% ============================================================
% 6. PI DE CORRIENTE EJE D
% ============================================================

fbw_d = PI_dq_params.fbw_d;                  % [Hz]
wd = 2*pi*fbw_d;

% Diseño usado habitualmente:
%
% Kp_d = wd*Ld
% Ki_d = wd*R

Kp_d = PI_dq_params.Kp_d_salient;
Ki_d = PI_dq_params.Ki_d_salient;

C_d = Kp_d + Ki_d/s;


%% ============================================================
% 7. PI DE CORRIENTE EJE Q
% ============================================================

fbw_q = PI_dq_params.fbw_q;                  % [Hz]
wq = 2*pi*fbw_q;

Kp_q = PI_dq_params.Kp_q_salient;
Ki_q = PI_dq_params.Ki_q_salient;

C_q = Kp_q + Ki_q/s;


%% ============================================================
% 8. PLANTAS ELECTRICAS
% ============================================================

% Aproximacion desacoplada:
%
% Id/Vd = 1/(Ld*s + R)
% Iq/Vq = 1/(Lq*s + R)

P_d = 1/(Ld*s + R);
P_q = 1/(Lq*s + R);


%% ============================================================
% 9. LAZO ABIERTO DE CORRIENTE
% ============================================================

L_d = C_d*P_d;
L_q = C_q*P_q;


%% ============================================================
% 10. LAZO CERRADO DE CORRIENTE
% ============================================================

% Id / Id_ref
T_d = L_d/(1 + L_d);

% Iq / Iq_ref
T_q = L_q/(1 + L_q);


%% ============================================================
% 11. PLANTA MECANICA
% ============================================================

% Wm / Te
%
% Gm(s) = 1/(J*s+B)

G_m = 1/(J*s + B);


%% ============================================================
% 12. MTPA
% ============================================================

% El MTPA es algebraico, por lo tanto no introduce fase dinamica.
%
% En primera aproximacion:
%
% Te_ref -> Te
%
% ganancia incremental = 1
%
% Si despues quieres calcular el Jacobiano exacto:
%
% dId_ref/dTe_ref
% dIq_ref/dTe_ref
%
% lo agregamos aqui.

K_MTPA = 1;

G_MTPA = K_MTPA;


%% ============================================================
% 13. MTPA + LAZO DE CORRIENTE
% ============================================================

% Si Id e Iq tienen BW parecidos, podemos tomar Tq como
% aproximacion de Te_ref -> Te.
%
% Para un modelo mas exacto se deben sumar las contribuciones
% de Id e Iq al torque.

G_Torque = G_MTPA*T_q;


%% ============================================================
% 14. CADENA DESDE Te_ref HASTA Wm
% ============================================================

G_Te_Wm = G_Torque*G_m;


%% ============================================================
% 15. CADENA DESDE Te_ref HASTA Wm_filtrada
% ============================================================

G_Te_WmFilt = G_Torque*G_m*H_Wm;


%% ============================================================
% 16. LAZO ABIERTO COMPLETO DE VELOCIDAD
% ============================================================

% error_Wm -> PI_Wm -> Te_ref -> ... -> Wm_filtrada

L_speed = C_w*G_Te_WmFilt;


%% ============================================================
% 17. LAZO CERRADO DE VELOCIDAD
% ============================================================

% Wm_ref -> Wm_filtrada

T_speed = L_speed/(1 + L_speed);


%% ============================================================
% 18. SENSIBILIDAD
% ============================================================

% Muy importante para analizar perturbaciones como cogging:
%
% S = 1/(1 + L)

S_speed = 1/(1 + L_speed);


%% ============================================================
% 19. RESPUESTA A UNA PERTURBACION DE TORQUE
% ============================================================

% Si el cogging entra mecanicamente como torque perturbador:
%
% Tcog -> Wm_filt
%
% signo negativo porque el cogging se opone al torque del motor.

G_cog_WmFilt_open = -G_m*H_Wm;

% Con PI de velocidad cerrando el lazo:
G_cog_WmFilt_closed = ...
    (-G_m*H_Wm)/(1 + L_speed);


%% ============================================================
% 20. FUNCION AUXILIAR PARA MOSTRAR RESULTADOS
% ============================================================

printFreqResponse('PI velocidad',C_w);

printFreqResponse('Filtro Wm',H_Wm);

printFreqResponse('PI Id',C_d);
printFreqResponse('Planta electrica d',P_d);
printFreqResponse('Lazo abierto Id',L_d);
printFreqResponse('Lazo cerrado Id',T_d);

printFreqResponse('PI Iq',C_q);
printFreqResponse('Planta electrica q',P_q);
printFreqResponse('Lazo abierto Iq',L_q);
printFreqResponse('Lazo cerrado Iq',T_q);

printFreqResponse('MTPA',G_MTPA);

printFreqResponse('MTPA + corriente',G_Torque);

printFreqResponse('Planta mecanica',G_m);

printFreqResponse('Te_ref -> Wm',G_Te_Wm);

printFreqResponse('Te_ref -> Wm filtrada',G_Te_WmFilt);

printFreqResponse('Lazo abierto velocidad',L_speed);

printFreqResponse('Lazo cerrado velocidad',T_speed);

printFreqResponse('Sensibilidad',S_speed);

printFreqResponse( ...
    'Cogging -> Wm filtrada (sin PI velocidad)', ...
    G_cog_WmFilt_open);

printFreqResponse( ...
    'Cogging -> Wm filtrada (lazo cerrado)', ...
    G_cog_WmFilt_closed);


%% ============================================================
% 21. TABLA RESUMEN
% ============================================================

names = {
    'PI Wm'
    'Filtro Wm'
    'PI Id'
    'Planta Id'
    'Closed-loop Id'
    'PI Iq'
    'Planta Iq'
    'Closed-loop Iq'
    'MTPA'
    'MTPA + corriente'
    'Mecanica'
    'TeRef -> Wm'
    'TeRef -> WmFilt'
    'Lazo velocidad'
    'Closed-loop velocidad'
    'Sensibilidad'
    'Cogging -> WmFilt OPEN'
    'Cogging -> WmFilt CLOSED'
    };

G_all = [
    C_w
    H_Wm
    C_d
    P_d
    T_d
    C_q
    P_q
    T_q
    G_MTPA
    G_Torque
    G_m
    G_Te_Wm
    G_Te_WmFilt
    L_speed
    T_speed
    S_speed
    G_cog_WmFilt_open
    G_cog_WmFilt_closed
    ];

mag_all = abs(G_all);

db_all = 20*log10(mag_all);

phase_all = angle(G_all)*180/pi;

ResultTable = table( ...
    names, ...
    mag_all, ...
    db_all, ...
    phase_all, ...
    'VariableNames', ...
    {'Bloque','Ganancia','Ganancia_dB','Fase_deg'});

disp(' ');
disp('============================================================');
fprintf(' RESULTADOS A %.3f Hz\n',f_cog);
fprintf(' omega = %.3f rad/s\n',w_cog);
disp('============================================================');

disp(ResultTable);


%% ============================================================
% 22. FASE ACUMULADA MANUAL
% ============================================================

phase_current = angle(G_Torque)*180/pi;
phase_mech    = angle(G_m)*180/pi;
phase_filter  = angle(H_Wm)*180/pi;

phase_total_unwrapped = ...
    phase_current + ...
    phase_mech + ...
    phase_filter;

fprintf('\n');
fprintf('============================================\n');
fprintf('FASE ACUMULADA Te_ref -> Wm_filt\n');
fprintf('============================================\n');

fprintf('MTPA + corriente : %9.3f deg\n',phase_current);
fprintf('Mecanica          : %9.3f deg\n',phase_mech);
fprintf('Filtro Wm         : %9.3f deg\n',phase_filter);

fprintf('--------------------------------------------\n');

fprintf('TOTAL SIN WRAP    : %9.3f deg\n', ...
    phase_total_unwrapped);

fprintf('TOTAL CON WRAP    : %9.3f deg\n', ...
    angle(G_Te_WmFilt)*180/pi);


%% ============================================================
% 23. FASE DEL PI DE VELOCIDAD
% ============================================================

phase_PI_w = angle(C_w)*180/pi;

fprintf('\n');
fprintf('PI velocidad      : %9.3f deg\n',phase_PI_w);

fprintf('Lazo velocidad    : %9.3f deg\n', ...
    angle(L_speed)*180/pi);


%% ============================================================
% FUNCION LOCAL
% ============================================================

function printFreqResponse(name,G)

    mag = abs(G);
    db = 20*log10(mag);
    phase = angle(G)*180/pi;

    fprintf('\n%-42s\n',name);
    fprintf('  Ganancia     = %12.6g\n',mag);
    fprintf('  Ganancia dB  = %12.4f dB\n',db);
    fprintf('  Fase         = %12.4f deg\n',phase);

end