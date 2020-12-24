clear; clc;
% system parameters
A = [-3, -2, 1; 1,2,1; 1, -1, -1];
B1 = [1 0.2 0]';
B2 = [0.4 1 0.6]';
C1 = [1.2 0.5 0.3];
D12 = 1;
% declare the variable
P = sdpvar(3);
Z = sdpvar(1);
F = sdpvar(1,3);
mu = sdpvar(1);
% describe the LMI
c1 = [A*P+P*A'+B2*F+F'*B2', P*C1'+F'*D12'; (P*C1'+F'*D12')', -eye(1)];
mat2 = [Z, B1';B1, P];
% declare the constraints
Fd = [c1<=0; mat2>=0; P>=0; Z>=0 ;trace(Z)<=mu];
optimize(Fd, mu);
% result
K = value(F)*inv(value(P))
H2_norm = sqrt(value(mu))