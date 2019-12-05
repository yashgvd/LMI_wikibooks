%% Declaring variables
X = sdpvar(3);
W = sdpvar(3,2);
Q=sdpvar(2);
rho = sdpvar(1);

%% DATA
A=[-3.2216 2.8762 0.6450;
-7.0385 -3.0859 -2.1610;
9.1159 4.5239 -6.3258];

B1=[2.6106 2.1274;
0.8399 2.0541;
0.2641 -1.8041];

B2=[0.1198;-0.0304;-0.0101];

C1=[3.0453 2.1022 1.9242;
4.4915 1.7021 -0.6568];

C2 =[1.9683 -1.4475 -2.2344;
-0.1943 0.6912 -1.3730];

D1=[-1.6956 1.7291;
0.5102 1.1855];

D2 =[0.0677
-0.0154];

eta=0.0001;
F=[];

%% Constraints
F1=[X*A+W*C1+(X*A+W*C1)' X*B2+W*D2;
(X*B2 + W*D2)' -1];

F2 = [-Q C2;C2' -X];

F=[F,X>0.0001*eye(3)];
F=[F,Q>0.0001*eye(2)];

F=[F,F1<eta*eye(size(F1))];
F=[F,F2<eta*eye(size(F2))];

F=[F,trace(Q)<rho];

%% Optimizing
optimize(F,rho);

%% Finding Observer Matrix
value(rho)
gamma = sqrt(value(rho))
L = inv(value(X))*value(W)