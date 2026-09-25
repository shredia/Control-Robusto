%% ========================================================================
%  ANÁLISIS CUANTITATIVO DE OBSERVABILIDAD
%  Motor paso a paso - modelo dq con saliencia
%
%  Estados:
%       x = [Id Iq Wm Theta_m]'
%
%  Salidas:
%       y = [Id Iq]'
%
%  Se calcula:
%
%   1) Jacobiano local A
%   2) Discretización Ad
%   3) Matriz de observabilidad de horizonte finito
%   4) Rango
%   5) Valores singulares
%   6) sigma_min
%   7) sigma_min/sigma_max
%   8) dirección débil v_min
%   9) Gramiano de observabilidad
%  10) autovalores del Gramiano
%  11) barrido respecto de velocidad
%  12) comparación con/sin saliencia
%
%  IMPORTANTE:
%
%  Las métricas cuantitativas se calculan sobre estados ESCALADOS.
%
% ========================================================================

clear;
clc;
close all;


%% ========================================================================
%  1. VARIABLES SIMBÓLICAS
% ========================================================================

syms Rs Ld Lq Ke Kt Jm Bm Nr real

syms Vd Vq real
syms Id Iq Wm Theta_m real

% Ángulo eléctrico utilizado por Park
syms Theta_hat_e real


%% ========================================================================
%  2. ÁNGULOS
% ========================================================================

Theta_e = Nr*Theta_m;

% Error angular
delta = Theta_e - Theta_hat_e;


%% ========================================================================
%  3. SALIENCIA
% ========================================================================

L0 = (Ld + Lq)/2;

L2 = (Ld - Lq)/2;


%% Matriz de inductancias en marco dq estimado

L11 = L0 + L2*cos(2*delta);

L22 = L0 - L2*cos(2*delta);

L12 = L2*sin(2*delta);


L_dq = [
    L11   L12
    L12   L22
];


%% ========================================================================
%  4. VECTORES ELÉCTRICOS
% ========================================================================

i_dq = [
    Id
    Iq
];


v_dq = [
    Vd
    Vq
];


%% ========================================================================
%  5. FLUJO PM
% ========================================================================
%
% NOTA:
% Se conserva aquí la misma convención del modelo que venimos
% analizando. Antes de utilizar resultados cuantitativos definitivos
% en la tesis debe verificarse la definición exacta de Ke respecto
% de la planta Simulink.
% ========================================================================

lambda_PM = Ke * [
    cos(delta)
    sin(delta)
];


%% ========================================================================
%  6. DERIVADAS ANGULARES
% ========================================================================

dL_dtheta = diff(L_dq,Theta_m);

dL_dt = dL_dtheta*Wm;


%% FEM

e_PM = jacobian(lambda_PM,Theta_m)*Wm;


%% ========================================================================
%  7. DINÁMICA ELÉCTRICA
% ========================================================================

rhs = ...
      v_dq ...
    - Rs*i_dq ...
    - dL_dt*i_dq ...
    - e_PM;


di_dq = L_dq \ rhs;


f1 = di_dq(1);

f2 = di_dq(2);


%% ========================================================================
%  8. TORQUE
% ========================================================================

T_PM = Kt*Iq;

T_rel = (Ld-Lq)*Id*Iq;

Te = T_PM + T_rel;


%% ========================================================================
%  9. DINÁMICA MECÁNICA
% ========================================================================

f3 = (Te-Bm*Wm)/Jm;

f4 = Wm;


%% ========================================================================
%  10. SISTEMA COMPLETO
% ========================================================================

x = [
    Id
    Iq
    Wm
    Theta_m
];


f = [
    f1
    f2
    f3
    f4
];


n = length(x);


%% ========================================================================
%  11. MATRIZ DE MEDICIÓN
% ========================================================================

C = [
    1 0 0 0
    0 1 0 0
];


%% ========================================================================
%  12. JACOBIANO
% ========================================================================

disp('Calculando Jacobiano simbólico...')

A_sym = jacobian(f,x);

disp('Jacobiano calculado.')


%% ========================================================================
%  13. PARÁMETROS
% ========================================================================
%
% Reemplaza posteriormente estos valores por los definitivos.
% ========================================================================

p.Rs = 2.5;

p.Ld = 0.0032;
p.Lq = 0.0081;

p.Ke = 0.045;
p.Kt = 0.045;

p.Jm = 0.01;
p.Bm = 0.11;

p.Nr = 50;


%% ========================================================================
%  14. PUNTO DE OPERACIÓN ELÉCTRICO
% ========================================================================

Id_op = 0;

Iq_op = 1;

Vd_op = 0;

Vq_op = 10;


%% ========================================================================
%  15. POSICIÓN
% ========================================================================

Theta_m_op = 0;

% Pequeño error angular inicial
Theta_hat_op = 0.05;


%% ========================================================================
%  16. TIEMPO DE MUESTREO
% ========================================================================
%
% Usa el tiempo de muestreo correspondiente al análisis del observador.
%
% Si tu EKF mecánico trabaja a 1 kHz:
%
%       Ts = 1e-3
%
% Si después analizas HFI/corrientes a mayor frecuencia,
% habrá que utilizar el Ts correspondiente.
% ========================================================================

Ts = 1e-3;


%% ========================================================================
%  17. HORIZONTE DE OBSERVACIÓN
% ========================================================================
%
% No estamos obligados a usar solamente N=n.
%
% Por ejemplo:
%
%       Nobs = 20
%
% con Ts = 1 ms equivale a observar 20 ms.
%
% Más adelante conviene estudiar sensibilidad a este horizonte.
% ========================================================================

Nobs = 20;


%% ========================================================================
%  18. ESCALAS DE LOS ESTADOS
% ========================================================================
%
% DEFINIMOS:
%
%       x = Sx*z
%
% donde z son estados normalizados.
%
% Las escalas deben representar magnitudes físicamente razonables.
%
% IMPORTANTE:
% Estos valores NO son universales.
% Debes justificarlos posteriormente.
% ========================================================================

I_scale = 1.0;       % [A]

W_scale = 10.0;      % [rad/s]

Theta_scale = 1.0;   % [rad]


Sx = diag([
    I_scale
    I_scale
    W_scale
    Theta_scale
]);


%% ========================================================================
%  19. VARIABLES PARA SUSTITUCIÓN
% ========================================================================

variables = [
    Rs
    Ld
    Lq
    Ke
    Kt
    Jm
    Bm
    Nr
    Vd
    Vq
    Id
    Iq
    Wm
    Theta_m
    Theta_hat_e
];


%% ========================================================================
%  20. FUNCIÓN PARA EVALUAR A
% ========================================================================
%
% Convertimos el Jacobiano simbólico en una función numérica.
%
% Esto hará que el barrido posterior sea mucho más rápido.
% ========================================================================

disp('Creando función numérica del Jacobiano...')

A_fun = matlabFunction( ...
    A_sym, ...
    'Vars', ...
    {Rs,Ld,Lq,Ke,Kt,Jm,Bm,Nr, ...
     Vd,Vq,Id,Iq,Wm,Theta_m,Theta_hat_e} ...
);

disp('Función creada.')


%% ========================================================================
%  21. PUNTO DE OPERACIÓN INICIAL
% ========================================================================

Wm_op = 0;


A_num = A_fun( ...
    p.Rs, ...
    p.Ld, ...
    p.Lq, ...
    p.Ke, ...
    p.Kt, ...
    p.Jm, ...
    p.Bm, ...
    p.Nr, ...
    Vd_op, ...
    Vq_op, ...
    Id_op, ...
    Iq_op, ...
    Wm_op, ...
    Theta_m_op, ...
    Theta_hat_op ...
);


%% ========================================================================
%  22. DISCRETIZACIÓN
% ========================================================================
%
% Sistema continuo:
%
%       dx/dt = A*x
%
% Sistema discreto:
%
%       x[k+1] = Ad*x[k]
%
% con:
%
%       Ad = exp(A*Ts)
%
% ========================================================================

Ad = expm(A_num*Ts);


%% ========================================================================
%  23. ESCALAMIENTO
% ========================================================================
%
% Si:
%
%       x = Sx*z
%
% entonces:
%
%       z[k+1] = inv(Sx)*Ad*Sx*z[k]
%
% y:
%
%       y = C*Sx*z
%
% ========================================================================

Ad_s = Sx \ (Ad*Sx);

C_s = C*Sx;


%% ========================================================================
%  24. MATRIZ DE OBSERVABILIDAD FINITA
% ========================================================================

O = [];

for k = 0:Nobs-1

    O = [
        O
        C_s*(Ad_s^k)
    ];

end


%% ========================================================================
%  25. RANGO
% ========================================================================

rank_O = rank(O);


%% ========================================================================
%  26. SVD
% ========================================================================

[U,S,V] = svd(O,'econ');

sv = diag(S);


sigma_max = max(sv);

sigma_min = min(sv);


%% ========================================================================
%  27. ÍNDICE NORMALIZADO
% ========================================================================

if sigma_max > 0

    eta = sigma_min/sigma_max;

else

    eta = 0;

end


%% ========================================================================
%  28. DIRECCIÓN DÉBIL
% ========================================================================

v_min = V(:,end);


%% ========================================================================
%  29. GRAMIANO DE OBSERVABILIDAD
% ========================================================================
%
% Para horizonte finito:
%
%       Wo = sum (Ad^k)' C'C Ad^k
%
% equivalentemente:
%
%       Wo = O'*O
%
% ========================================================================

Wo = O'*O;


%% ========================================================================
%  30. AUTOVALORES DEL GRAMIANO
% ========================================================================

lambda_W = eig(Wo);

lambda_W = sort(real(lambda_W),'descend');


lambda_W_max = max(lambda_W);

lambda_W_min = min(lambda_W);


%% ========================================================================
%  31. RESULTADOS PUNTO DE OPERACIÓN
% ========================================================================

fprintf('\n============================================\n')
fprintf('OBSERVABILIDAD - PUNTO DE OPERACIÓN\n')
fprintf('============================================\n')

fprintf('Wm                = %.4f rad/s\n',Wm_op);

fprintf('Rank(O)           = %d / %d\n', ...
    rank_O,n);

fprintf('Sigma max         = %.6e\n', ...
    sigma_max);

fprintf('Sigma min         = %.6e\n', ...
    sigma_min);

fprintf('Sigma min/max     = %.6e\n', ...
    eta);

fprintf('Lambda min(Wo)    = %.6e\n', ...
    lambda_W_min);

fprintf('Lambda max(Wo)    = %.6e\n', ...
    lambda_W_max);


%% ========================================================================
%  32. COMPROBACIÓN RELACIÓN GRAMIANO-SVD
% ========================================================================
%
% Teóricamente:
%
%       lambda_i(Wo) = sigma_i(O)^2
%
% ========================================================================

fprintf('\nComprobación:\n');

fprintf('sigma_min^2       = %.6e\n', ...
    sigma_min^2);

fprintf('lambda_min(Wo)    = %.6e\n', ...
    lambda_W_min);


%% ========================================================================
%  33. DIRECCIÓN MENOS OBSERVABLE
% ========================================================================

state_names = {
    'Id'
    'Iq'
    'Wm'
    'Theta_m'
};


fprintf('\nDirección menos observable:\n');


for i = 1:n

    fprintf('%-10s : %+12.6e\n', ...
        state_names{i}, ...
        v_min(i));

end


%% ========================================================================
%  34. PARTICIPACIÓN NORMALIZADA
% ========================================================================

participation = abs(v_min);

participation = participation/max(participation);


fprintf('\nParticipación dirección débil:\n');


for i = 1:n

    fprintf('%-10s : %.4f\n', ...
        state_names{i}, ...
        participation(i));

end


%% ========================================================================
%  35. BARRIDO DE VELOCIDAD
% ========================================================================
%
% Vamos desde velocidad cero hasta 10 rad/s.
%
% Incluimos mayor densidad cerca de cero.
% ========================================================================

W_vector = [
    0 ...
    logspace(-4,1,150)
];


Nw = length(W_vector);


%% Reservar memoria

rank_sal = zeros(Nw,1);

sigmaMin_sal = zeros(Nw,1);

eta_sal = zeros(Nw,1);

lambdaMin_sal = zeros(Nw,1);


rank_iso = zeros(Nw,1);

sigmaMin_iso = zeros(Nw,1);

eta_iso = zeros(Nw,1);

lambdaMin_iso = zeros(Nw,1);


%% ========================================================================
%  36. MOTOR ISOTRÓPICO
% ========================================================================

Liso = (p.Ld+p.Lq)/2;


%% ========================================================================
%  37. LOOP DE VELOCIDAD
% ========================================================================

disp(' ')
disp('Realizando barrido de velocidad...')


for j = 1:Nw

    W_test = W_vector(j);


    %% ================================================================
    % CASO A: MOTOR CON SALIENCIA
    % ================================================================

    A_test = A_fun( ...
        p.Rs, ...
        p.Ld, ...
        p.Lq, ...
        p.Ke, ...
        p.Kt, ...
        p.Jm, ...
        p.Bm, ...
        p.Nr, ...
        Vd_op, ...
        Vq_op, ...
        Id_op, ...
        Iq_op, ...
        W_test, ...
        Theta_m_op, ...
        Theta_hat_op ...
    );


    %% Discretización

    Ad_test = expm(A_test*Ts);


    %% Escalamiento

    Ad_test_s = Sx \ (Ad_test*Sx);


    %% Matriz O

    O_test = [];

    for k = 0:Nobs-1

        O_test = [
            O_test
            C_s*(Ad_test_s^k)
        ];

    end


    %% SVD

    s_test = svd(O_test);


    rank_sal(j) = rank(O_test);

    sigmaMin_sal(j) = min(s_test);

    eta_sal(j) = ...
        min(s_test)/max(s_test);


    %% Gramiano

    Wo_test = O_test'*O_test;

    eig_test = eig(Wo_test);

    lambdaMin_sal(j) = ...
        max(0,min(real(eig_test)));


    %% ================================================================
    % CASO B: MOTOR SIN SALIENCIA
    % ================================================================

    A_test_iso = A_fun( ...
        p.Rs, ...
        Liso, ...
        Liso, ...
        p.Ke, ...
        p.Kt, ...
        p.Jm, ...
        p.Bm, ...
        p.Nr, ...
        Vd_op, ...
        Vq_op, ...
        Id_op, ...
        Iq_op, ...
        W_test, ...
        Theta_m_op, ...
        Theta_hat_op ...
    );


    %% Discretización

    Ad_test_iso = expm(A_test_iso*Ts);


    %% Escalamiento

    Ad_test_iso_s = ...
        Sx \ (Ad_test_iso*Sx);


    %% Observabilidad

    O_test_iso = [];

    for k = 0:Nobs-1

        O_test_iso = [
            O_test_iso
            C_s*(Ad_test_iso_s^k)
        ];

    end


    %% SVD

    s_iso = svd(O_test_iso);


    rank_iso(j) = rank(O_test_iso);

    sigmaMin_iso(j) = min(s_iso);

    eta_iso(j) = ...
        min(s_iso)/max(s_iso);


    %% Gramiano

    Wo_iso = O_test_iso'*O_test_iso;

    eig_iso = eig(Wo_iso);

    lambdaMin_iso(j) = ...
        max(0,min(real(eig_iso)));

end


disp('Barrido terminado.')


%% ========================================================================
%  38. GRÁFICA: RANGO VS VELOCIDAD
% ========================================================================

figure

semilogx( ...
    max(W_vector,1e-4), ...
    rank_sal, ...
    'LineWidth',1.5)

hold on

semilogx( ...
    max(W_vector,1e-4), ...
    rank_iso, ...
    '--', ...
    'LineWidth',1.5)

grid on

xlabel('\omega_m [rad/s]')

ylabel('rank(\mathcal{O})')

title('Rango de observabilidad vs velocidad')

legend( ...
    'Con saliencia', ...
    'Sin saliencia', ...
    'Location','best')


%% ========================================================================
%  39. GRÁFICA: SIGMA MIN
% ========================================================================

figure

loglog( ...
    max(W_vector,1e-4), ...
    max(sigmaMin_sal,eps), ...
    'LineWidth',1.5)

hold on

loglog( ...
    max(W_vector,1e-4), ...
    max(sigmaMin_iso,eps), ...
    '--', ...
    'LineWidth',1.5)

grid on

xlabel('\omega_m [rad/s]')

ylabel('\sigma_{min}(\mathcal{O})')

title('Menor valor singular vs velocidad')

legend( ...
    'Con saliencia', ...
    'Sin saliencia', ...
    'Location','best')


%% ========================================================================
%  40. GRÁFICA: ÍNDICE SIGMA MIN / SIGMA MAX
% ========================================================================

figure

loglog( ...
    max(W_vector,1e-4), ...
    max(eta_sal,eps), ...
    'LineWidth',1.5)

hold on

loglog( ...
    max(W_vector,1e-4), ...
    max(eta_iso,eps), ...
    '--', ...
    'LineWidth',1.5)

grid on

xlabel('\omega_m [rad/s]')

ylabel('\sigma_{min}/\sigma_{max}')

title('Índice relativo de observabilidad')

legend( ...
    'Con saliencia', ...
    'Sin saliencia', ...
    'Location','best')


%% ========================================================================
%  41. GRÁFICA: LAMBDA MIN DEL GRAMIANO
% ========================================================================

figure

loglog( ...
    max(W_vector,1e-4), ...
    max(lambdaMin_sal,eps), ...
    'LineWidth',1.5)

hold on

loglog( ...
    max(W_vector,1e-4), ...
    max(lambdaMin_iso,eps), ...
    '--', ...
    'LineWidth',1.5)

grid on

xlabel('\omega_m [rad/s]')

ylabel('\lambda_{min}(W_o)')

title('Mínimo autovalor del Gramiano de observabilidad')

legend( ...
    'Con saliencia', ...
    'Sin saliencia', ...
    'Location','best')


%% ========================================================================
%  42. RESULTADO ESPECIAL EN VELOCIDAD CERO
% ========================================================================

fprintf('\n============================================\n')
fprintf('RESULTADOS EN VELOCIDAD CERO\n')
fprintf('============================================\n')


fprintf('\nCON SALIENCIA\n')

fprintf('Rank              = %d/%d\n', ...
    rank_sal(1),n);

fprintf('Sigma min         = %.6e\n', ...
    sigmaMin_sal(1));

fprintf('Sigma min/max     = %.6e\n', ...
    eta_sal(1));

fprintf('Lambda min Wo     = %.6e\n', ...
    lambdaMin_sal(1));


fprintf('\nSIN SALIENCIA\n')

fprintf('Rank              = %d/%d\n', ...
    rank_iso(1),n);

fprintf('Sigma min         = %.6e\n', ...
    sigmaMin_iso(1));

fprintf('Sigma min/max     = %.6e\n', ...
    eta_iso(1));

fprintf('Lambda min Wo     = %.6e\n', ...
    lambdaMin_iso(1));


%% ========================================================================
%  43. FIN
% ========================================================================

disp(' ')
disp('============================================')
disp('ANÁLISIS COMPLETO FINALIZADO')
disp('============================================')