clear all
clc 

%% Declaring variables
P = sdpvar(2,2);
K = sdpvar(1);
X = sdpvar(2,2);
%% DATA
A=[0 1; 1 0];

B=[1;0];


C=[1 10];

eta=0.0001;


%% Constraints
F = [X>(1e-6)*eye(2)];
F = [F;[A'*X+X*A-P*B*B'*X-X*B*B'*P+X*B*B'*X P*B+C'*K';(P*B+C'*K')' -1 ]<eta*eye(3)]; %The lyapunov equation
F = [F;P>(1e-6)*eye(2)];


optimize(F,P);

value(P)
value(K)