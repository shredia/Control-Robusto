clear;
clc;

%% Estados
syms ia ib th w real
x = [ia; ib; th; w];

%% Salida
h = [ia; ib];

%% DEFINIR f1..f4 COMO FUNCIONES DE TODOS LOS ESTADOS
syms f1(ia,ib,th,w) f2(ia,ib,th,w) f3(ia,ib,th,w) f4(ia,ib,th,w)
f = [f1; f2; f3; f4];

%% Operador de Lie
Lie = @(g) jacobian(g,x)*f;

%% Derivadas de Lie
Lf_h  = Lie(h);
Lf2_h = Lie(Lf_h);
Lf3_h = Lie(Lf2_h);

%% Matriz de Observabilidad
O1 = jacobian(h,x);
O2 = jacobian(Lf_h,x);
O3 = jacobian(Lf2_h,x);
O4 = jacobian(Lf3_h,x);

ObservabilityMatrix = [O1; O2; O3; O4];
OBS_simplified = simplify( subs(ObservabilityMatrix, {
    diff(f1,ib), diff(f1,th), diff(f2,ia), diff(f3,ib)
}, {
    0, 0, 0, 0
}))