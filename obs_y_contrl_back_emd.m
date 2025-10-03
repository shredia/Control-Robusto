clear;
clc;
R =  3.4; %% [Ohms]
L = 6/1000; %% [L]
Nr = 50; %%Relación entre torque electrico y mecánico
Km = 0.275; %%constante de torque xd
Tl = 1; %%Carga del motor
B = 5/10000; %%roce
Wm = 1; %vueltas por segundo
Iq = Tl/Km; %%Corriente Q
Id = 0; %%Corriente D
Vd = 0; 
Vq = Tl*R/Km;
J = 8/1000; %%inercia del motor
We = 5000;
Ia = 1.5;
Ib = 1.5;
Ke = 0.0055;
Kt = Nr*Ke;
theta = 45;


A = [(-R/L) 0 (Ke*sin(theta)/L) (Ke*We*cos(theta)/L) 0;
     0 (-R/L)  (-Ke*cos(theta)/L) ((Ke*We*sin(theta))/L) 0;
     0 0 0 1 0;
      (-Kt*sin(theta)/J) (Kt*cos(theta)/J) -B/J (Kt/J)*(-Ib*sin(theta)-Ia*cos(theta)) -1/J;
      0 0 0 0 0]

B = [1/L 0;
      0  1/L;
      0   0;
      0   0;
      0   0]

C = [1 0 0 0 0;0 1 0 0 0 ];


D = [0 0];

matriz_observabilidad = [C;C*A;C*A*A;C*A*A*A]
observavilidad = rank(matriz_observabilidad)

matriz_controlabilidad = [B A*B A*A*B A*A*A*B]
controlabilidad = rank(matriz_controlabilidad)

[U,S,V] = svd(matriz_observabilidad);
diag(S)   % te da los valores singulares

% Polos deseados del observador (ajusta si hace falta)


% Ganancia L continua (5x2)
eig(A)
L = place(A', C', 10*eig(A)).';
fprintf('L = [\n');
for i=1:size(L,1)
    fprintf('  ');
    fprintf('%12.6f ', L(i,:));
    fprintf(';\n');
end
fprintf('];\n');
