clear; clc;
I = 1;
Vdc = 12;
Kt = 0.33;
R = 3.4; 
L = 6/1000; % H
P = 50;

Aa = (L*I)^2 + (Kt/P)^2;
Ba = (-2*R*I*Kt/P);
Ca = (R*I)^2 - Vdc^2;

S1a = (-Ba + sqrt((Ba)^2 - 4*Aa*Ca)) / (2*Aa*2*pi)
S2a = (-Ba - sqrt((Ba)^2 - 4*Aa*Ca)) / (2*Aa*2*pi)

Ab = (L*I)^2 + (Kt/P)^2;
Bb = (2*R*I*Kt/P);
Cb = (R*I)^2 - Vdc^2;

S1b = (-Bb + sqrt((Bb)^2 - 4*Ab*Cb)) / (2*Ab*2*pi)
S2b = (-Bb - sqrt((Bb)^2 - 4*Ab*Cb)) / (2*Ab*2*pi)

%%S1a 471 hz S1b 381 hz A 24 volts.
%%S1a 255 hz S1b 165 hz A 24 volts.