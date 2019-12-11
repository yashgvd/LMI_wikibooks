clc;clear all
A=[-2 0 1;
    0 -3 0; 
    1 0 -2];
Ad=[-1 1 1;
    2 -1 1
    0 0 -1];

beta = sdpvar(1);
X=sdpvar(3);
db=0.1;

Phi_X=X*(A+Ad)'+(A+Ad)*X+db*Ad*Ad';

mat1=[Phi_X db*(A*X)' db*X*(Ad');
   db*(A*X) -db*beta*eye(3) zeros(3);
   db*Ad*X  zeros(3) -db*(1-beta)*eye(3)];

F=[X>=0.0001*eye(size(A,1))];
F=[F;;0<=beta<=1];
F=[F;mat1<=0];

optimize(F);

X=value(X)
beta=value(beta)