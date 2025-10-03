clear;
clc;
R =  3.4; %% [Ohms]
L = 6/1000; %% [L]
Nr = 50; %%Relación entre torque electrico y mecánico
Km = 10; %%que era esto? xd
Tl = 1; %%Carga del motor
B = 5/10000; %%roce
Wm = 1; %vueltas por segundo
Iq = Tl/Km; %%Corriente Q
Id = 0; %%Corriente D
Vd = 0; 
Vq = Tl*R/Km;
J = 8/1000; %%inercia del motor
We = 0;
tetha = 90;



A = [-R/L                     0                  (1/L)      0       0    0;
       0                    -R/L                   0     (-1/L)     0    0;
       0                      0                    We       0       0    0;
       0                      0                   -We       0       0    0
((Km*cos(tetha))/J)   -((Km*cos(tetha))/J)         0        0    (-B/J)  1;
        0                     0                    0        0       0    0];

B = [1/L 0;
      0  1/L;
      0   0;
      0   0;
      0   0;
      0   0];

C = [1 0 0 0 0 0;
     0 1 0 0 0 0];


D = [0 0];

matriz_observabilidad = [C;C*A;C*A*A;C*A*A*A]
observavilidad = rank(matriz_observabilidad)

matriz_controlabilidad = [B A*B A*A*B A*A*A*B]
controlabilidad = rank(matriz_controlabilidad)
