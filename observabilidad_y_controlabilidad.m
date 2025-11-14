%% === Sistema no lineal y cálculo de Jacobianos ===
clear; clc;
syms R L Kt Va Vb Ia Ib We Wm Tx J Th_m Nr B % Variables simbólicas

% --- Ecuaciones del sistema ---
f1 = (Va - R*Ia + Kt*Wm*sin(Nr*Th_m))/L;                % d(Id)/dt
f2 = (Vb - R*Ib - Kt*Wm*cos(Nr*Th_m))/L;        % d(Iq)/dt
f3 = (Kt*(Ib*cos(Nr*Th_m)-Ia*sin(Nr*Th_m)) -B*Wm - Tx)/J;                         % d(Wm)/dt
f4 = Wm;                                     % d(Th)/dt  (si se usa velocidad eléctrica)
f5 = 0;                                      % d(Tx)/dt (constante o perturbación lenta)
f = [f1; f2; f3; f4; f5];                    % Vector de funciones
x = [Ia; Ib; Wm; Th_m; Tx];                    % Vector de estados
u = [Va; Vb];                                % Vector de entradas

% --- Jacobianos ---
A = jacobian(f, x);  % Matriz A = ∂f/∂x
B = jacobian(f, u);  % Matriz B = ∂f/∂u

% --- Matriz de observación (por ejemplo, medición de Id e Iq) ---
C = [1 0 0 0 0;
     0 1 0 0 0];

D = zeros(2,2);

% --- Resultados ---
disp('Matriz A = ');
pretty(A)
disp('Matriz B = ');
pretty(B)
%disp('Matriz C = ');
%disp(C)
%disp('Matriz D = ');
%disp(D)

% --- Ver observabilidad y controlabilidad ---
matriz_observabilidad = [C; C*A; C*A^2; C*A^3;C*A^4];
%pretty(matriz_observabilidad);
obs_rank = rank(matriz_observabilidad);
disp(['Rango observabilidad: ', num2str(obs_rank)]);

matriz_controlabilidad = [B A*B A^2*B A^3*B A^4*B];
%pretty(matriz_controlabilidad);
ctrl_rank = rank(matriz_controlabilidad);
disp(['Rango controlabilidad: ', num2str(ctrl_rank)]);

