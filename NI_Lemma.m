clear all
clc

%% Data
A=[1 2; 3 4];
B=[1;6];
C=[4 7];
D=[5];
P=sdpvar(2);

%% Constraining and Optimizing
F1=[P*A+A'*P P*B-A'*C';(P*B-A'*C')' -(C*B+B'*C')]; 

F=[P>.001*eye(2)];
F=[F;F1<.001*eye(3)];

optimize(F,P)
value(P)