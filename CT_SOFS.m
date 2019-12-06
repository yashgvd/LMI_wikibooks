clear all
clc 

%% Declaring variables
P = sdpvar(3,3);
K = sdpvar(1,2);

%% DATA
A=[-4.701 1 0;-8.2986 0 0; 1 0 0];

B=[-0.0721;15.0218;0];


C=[1 0 0; 0 0 1];

eta=0.0001;

%% Constraints
F1=[A'*P+P*A-P*B*B'*P P*B+C'*K';
K*C+B'*P -1];

F=[P>(1e-6)*eye(3)];

F=[F;F1<0];

%% Optimizing
options=sdpsettings('solver','sedumi');
optimize(F,P);
%%
value(K)