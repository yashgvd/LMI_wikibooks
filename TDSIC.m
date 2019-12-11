clc;clear all
A=[-2 0 1;
    0 -3 0; 
    1 0 -2];
Ad=[-1 1 1;
    2 -1 1
    0 0 -1];
P=sdpvar(size(A,1));
S=sdpvar(size(A,1));
mat1=[A'*P+P*A+S   P*Ad;
   P*Ad' -S];
F=[P>.00001*eye(size(A,1))];
F=[F;mat1<0];
optimize(F);
value(P)