%% Discretización exacta (ZOH) del modelo iα–iβ–eα–eβ + chequeos
% Autor: tú :)
% Requiere MATLAB/Octave con expm. Probado con ambos.

clear; clc;

%% === Parámetros NUMÉRICOS de ejemplo (pon los tuyos) ===
R   = 3.4;           % Ohm
L   = 0.01;          % H
We  = 0.532;         % rad/s   (ejemplo nominal)  <<--- PON TU VALOR
Ts  = 1.0;           % s       (periodo de muestreo) <<--- PON TU VALOR

%% === Construcción de A, B, C (continuo) ===
% x = [ i_alpha ; i_beta ; e_alpha ; e_beta ]
% u = [ v_alpha ; v_beta ]
A = [ -R/L,    0,   -1/L,   0;
        0,   -R/L,    0,   -1/L;
        0,     0,     0,   We;
        0,     0,   -We,    0 ]

B = [ 1/L, 0;
      0,  1/L;
      0,    0;
      0,    0 ]

C = [ 1, 0, 0, 0;
      0, 1, 0, 0 ]

n = size(A,1);
m = size(B,2);

%% === A) Discretización exacta por expm(A*Ts) + Bd con identidad ZOH ===
Ad_expm = expm(A*Ts);

% Bd = ∫_0^Ts e^{A τ} B dτ = A^{-1} (Ad - I) B    (ZOH)
% Usa mldivide por estabilidad: A \ ((Ad - I)*B)
Bd_expm = A \ ((Ad_expm - eye(n)) * B);

