%% === Sistema no lineal y cálculo de Jacobianos ===
clear; clc;
syms R L Kt Va Vb Ia Ib We Wm Tl J Th_m Nr B a b Wa Wb Ke Tx Te% Variables simbólicas

% --- Ecuaciones del sistema ---
f1_1 = (Va - R*Ia + Ke*Wm*sin(Nr*Th_m))/L;                % d(Id)/dt
f2_1 = (Vb - R*Ib - Ke*Wm*cos(Nr*Th_m))/L;        % d(Iq)/dt
f3_1 = (Kt*(Ib*cos(Nr*Th_m)-Ia*sin(Nr*Th_m)) -Tx -B*Wm)/J;                         % d(Wm)/dt
f4_1 = Wm; 

f1_2 = (Te-Tx)/J;                % d(Id)/dt
f2_2 = (0);        % d(Iq)/dt



% d(Th)/dt  (si se usa velocidad eléctrica)                                     % d(Tx)/dt (constante o perturbación lenta)
f_1 = [f1_1; f2_1; f3_1; f4_1];                    % Vector de funciones
x_1 = [Ia; Ib; Wm; Th_m];                    % Vector de estados
u_1 = [Va; Vb];                                % Vector de entradas
% d(Th)/dt  (si se usa velocidad eléctrica)                                     % d(Tx)/dt (constante o perturbación lenta)
f_2 = [f1_2; f2_2;];                    % Vector de funciones
x_2 = [Wm; Tx];                    % Vector de estados
u_2 = Te;                                % Vector de entradas

% --- Jacobianos ---
A_1 = jacobian(f_1,x_1);  % Matriz A = ∂f/∂x
B_1 = jacobian(f_1,u_1);  % Matriz B = ∂f/∂u

A_2 = jacobian(f_2,x_2);  % Matriz A = ∂f/∂x
B_2 = jacobian(f_2, u_2);  % Matriz B = ∂f/∂u

% --- Matriz de observación (por ejemplo, medición de Id e Iq) ---
% --- Matriz de observación (por ejemplo, medición de Id e Iq) ---
C_1 = [1 0 0 0  ;
       0 1 0 0 ];

C_2 = [1 0 ];

D_1 = zeros(2,2);
D_2 = zeros(1,1);

% --- Resultados ---
disp('Matriz A_1 = ');
pretty(A_1)
disp('Matriz B_1 = ');
pretty(B_1)
disp('Matriz A_2 = ');
pretty(A_2)
disp('Matriz B_2 = ');
pretty(B_2)
%disp('Matriz C = ');
%disp(C)
%disp('Matriz D = ');
%disp(D)

% --- Ver observabilidad y controlabilidad ---
matriz_observabilidad_1 = [C_1; C_1*A_1; C_1*A_1^2; C_1*A_1^3]
%pretty(matriz_observabilidad);

obs_rank_1 = rank(matriz_observabilidad_1);
disp(['Rango observabilidad: ', num2str(obs_rank_1)]);

matriz_controlabilidad_1 = [B_1 A_1*B_1 A_1^2*B_1 A_1^3*B_1 A_1^4*B_1];
%pretty(matriz_controlabilidad);
ctrl_rank_1 = rank(matriz_controlabilidad_1);
disp(['Rango controlabilidad: ', num2str(ctrl_rank_1)]);
v_no_obs_1 = null(matriz_observabilidad_1)   % devuelve una base del subespacio no observable

% --- Ver observabilidad y controlabilidad ---
matriz_observabilidad_2 = [C_2; C_2*A_2;]
%pretty(matriz_observabilidad);

obs_rank_2 = rank(matriz_observabilidad_2);
disp(['Rango observabilidad: ', num2str(obs_rank_2)]);

matriz_controlabilidad_2 = [B_2 A_2*B_2 ];
%pretty(matriz_controlabilidad);
ctrl_rank_2 = rank(matriz_controlabilidad_2);
disp(['Rango controlabilidad: ', num2str(ctrl_rank_2)]);
v_no_obs_2 = null(matriz_observabilidad_2)   % devuelve una base del subespacio no observable

