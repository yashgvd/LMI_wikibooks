clear all
clc 

%% Declaring variables
P = sdpvar(4,4);
K = sdpvar(2,2);

%% DATA
A=[0.5 0 0.2 1.0;
0 0.3 0 0.1;
0.01 0.1 -0.5  0;
0.1 0 -0.1 -1.0];

B=[1 0; 0 1; 0 0 ; -1 0];


C=[1 0 0 1;1 0 1 1];

eta=0.0001;

%% Constraints
F1=[-P (A+B*K*C)*P;
((A+B*K*C)*P)' -P];

F=[P>(1e-6)*eye(4)];

F=[F;F1<0];

%% Optimizing
options=sdpsettings('solver','sedumi');
optimize(F,P);
%%
value(K)