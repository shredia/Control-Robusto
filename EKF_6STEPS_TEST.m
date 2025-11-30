clear;
clc;

% ====== VECTORES DE ESTADO ======
syms x1 x2 x3 x4 x5 real
x_hat_less     = [x1; x2; x3; x4; x5];
x_dot_hat_less = sym('xdot', [5 1], 'real');   % xdot1...xdot5

% ====== MATRICES ======
sum_less = sym('sum_less', [5 5], 'real');
sum_more = sym('sum_more', [5 5], 'real');
sum_v = sym('sum_v', [2 2], 'real');


% ====== VECTOR DE SALIDA ======
y_hat = sym('yhat', [2 1], 'real');            % yhat1, yhat2

% ====== Mostrar resultados ======
disp('x_hat_less ='); disp(x_hat_less);
disp('x_dot_hat_less ='); disp(x_dot_hat_less);
disp('y_hat ='); disp(y_hat);
C = [1 0 0 0 0; 0 1 0 0 0];
CT = C';
Lk = simplify(sum_less*CT*inv(C*sum_less*CT + sum_v));
disp('lK ='); disp(simplify(Lk));
C*sum_less*CT
