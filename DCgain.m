clear all
clc

%% Data
A=[1 2; 3 4];
B=[1;6];
C=[4 7];
D=[5];
gamma=sdpvar(1);

%% Constraining and Optimizing
F1=[gamma -C*inv(A)*B+D;(-C*inv(A)*B+D)' gamma]; 


F=[F1>.001*eye(2)];

optimize(F,gamma)
value(gamma)
%% DC Gain
G=D-C*inv(A)*B