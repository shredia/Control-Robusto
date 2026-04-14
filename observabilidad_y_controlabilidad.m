%% === Sistema no lineal y cálculo de Jacobianos ===
clear; clc;

 syms Rs R Ld Lq Kt Vd Vq Id Iq We Wm Tl J Th_m Nr B a b Wa Wb Ke Tx Te L0 L2 La Lb% Variables simbólicas

f1_1 = (Vd - Rs*Id + Nr*Wm*Lq*Iq)/Ld;                % d(Id)/dt
f2_1 = (Vq - Rs*Iq - Nr*Wm*Ld*Id-Ke*Wm)/Lq;        % d(Iq)/dt
f3_1 = (Kt*(Iq) + (Ld-Lq)*Id*Iq -B*Wm)/J;                         % d(Wm)/dt
f4_1 = Wm; 
f5_1 = 0;
% % d(Th)/dt  (si se usa velocidad eléctrica)                                     % d(Tx)/dt (constante o perturbación lenta)
 f_1 = [f1_1; f2_1; f3_1; f4_1];                    % Vector de funciones
 x_1 = [Id; Iq; Wm; Th_m];                    % Vector de estados
 u_1 = [Vd; Vq];                                % Vector de entradas
 
% % --- Jacobianos ---
 A_1 = jacobian(f_1,x_1);  % Matriz A = ∂f/∂x
 B_1 = jacobian(f_1,u_1);  % Matriz B = ∂f/∂u

% % --- Matriz de observación (por ejemplo, medición de Id e Iq) ---
C_1 = [1 0 0 0  ;
       0 1 0 0 ;
       0 0 0 1 ];


D_1 = zeros(2,2);

% % --- Resultados ---
disp('Matriz A_1 = ');
pretty(A_1)
disp('Matriz B_1 = ');
pretty(B_1)

% 
% --- Ver observabilidad y controlabilidad ---
% --- Construcción de la Matriz de Observabilidad ---
O = [C_1; 
     C_1*A_1; 
     C_1*A_1^2];

% --- SUSTITUCIÓN NUMÉRICA (Crítico para que rank y cond funcionen) ---
% Definimos un struct con tus valores medidos
params.Rs = 2.5; 
params.Ld = 0.0032; 
params.Lq = 0.0081; 
params.Ke = 0.045; 
params.Kt = 0.045; 
params.J  = 0.01; 
params.B  = 0.11; 
params.Nr = 50; 
params.Wm = 10; % Tu velocidad de interés (0.1 rad/s o RPM)
params.Id = 0; 
params.Iq = 1; 
params.Tx = 0.2;
params.Vd = 0;
params.Vq = 10;

% Convertimos la matriz simbólica a numérica
O_num = double(subs(O, fieldnames(params), struct2cell(params)));

% --- Resultados Finales ---
obs_rank = rank(O_num);
cond_num = cond(O_num);

fprintf('\n--- RESULTADOS DEL ANÁLISIS ---\n');
fprintf('Rango de Observabilidad: %d\n', obs_rank);
fprintf('Número de Condición: %.2e\n', cond_num);

if obs_rank < 5
    disp('El sistema NO es completamente observable.');
    v_no_obs = null(O_num);
    disp('Base del subespacio no observable:');
    disp(v_no_obs);
else
    disp('El sistema es COMPLETAMENTE OBSERVABLE.');
end

s = svd(O_num);
disp('Valores singulares de la matriz de observabilidad:');
disp(s);










% syms R_a R_b L Kt Va Vb Ia Ib We Wm Tl J Th_m Nr B a b Wa Wb Ke Tx Te L0 L2 La Lb% Variables simbólicas
% La = L0 + L2*cos(2*Th_m*Nr);
% Lb = L0 - L2*cos(2*Th_m*Nr);
% % --- Ecuaciones del sistema ---
% f1_1 = (Va - R_a*Ia + Ke*Wm*sin(Nr*Th_m))/La;                % d(Id)/dt
% f2_1 = (Vb - R_b*Ib - Ke*Wm*cos(Nr*Th_m))/Lb;        % d(Iq)/dt
% f3_1 = (Kt*(Ib*cos(Nr*Th_m)-Ia*sin(Nr*Th_m)) -Tx -B*Wm)/J;                         % d(Wm)/dt
% f4_1 = Wm; 
% 
% f1_2 = (Te-Tx)/J;                % d(Id)/dt
% f2_2 = (0);        % d(Iq)/dt
% 
% 
% 
% % d(Th)/dt  (si se usa velocidad eléctrica)                                     % d(Tx)/dt (constante o perturbación lenta)
% f_1 = [f1_1; f2_1; f3_1; f4_1];                    % Vector de funciones
% x_1 = [Ia; Ib; Wm; Th_m];                    % Vector de estados
% u_1 = [Va; Vb];                                % Vector de entradas
% 
% % --- Jacobianos ---
% A_1 = jacobian(f_1,x_1);  % Matriz A = ∂f/∂x
% B_1 = jacobian(f_1,u_1);  % Matriz B = ∂f/∂u
% 
% 
% % --- Matriz de observación (por ejemplo, medición de Id e Iq) ---
% % --- Matriz de observación (por ejemplo, medición de Id e Iq) ---
% C_1 = [1 0 0 0  ;
%        0 1 0 0 ];
% 
% 
% 
% D_1 = zeros(2,2);
% 
% 
% % --- Resultados ---
% disp('Matriz A_1 = ');
% pretty(A_1)
% disp('Matriz B_1 = ');
% pretty(B_1)
% 
% 
% % --- Ver observabilidad y controlabilidad ---
% matriz_observabilidad_1 = [C_1; C_1*A_1; C_1*A_1^2; C_1*A_1^3]
% %pretty(matriz_observabilidad);
% 
% obs_rank_1 = rank(matriz_observabilidad_1);
% disp(['Rango observabilidad: ', num2str(obs_rank_1)]);
% 
% matriz_controlabilidad_1 = [B_1 A_1*B_1 A_1^2*B_1 A_1^3*B_1 A_1^4*B_1];
% %pretty(matriz_controlabilidad);
% ctrl_rank_1 = rank(matriz_controlabilidad_1);
% disp(['Rango controlabilidad: ', num2str(ctrl_rank_1)]);
% s
% v_no_obs_1 = null(matriz_observabilidad_1)   % devuelve una base del subespacio no observable
% 

