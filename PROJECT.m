m=1140;
Iz=1436.24;
Lr=1.165;
Lf=1.165;
Cf=155494.663;
Cr=155494.663;
Vx=17.3205;

A=[ 0 1 0 0;
    0 -2*(Cf+Cr)/(m*Vx) 2*(Cf+Cr)/m -Vx-2*(Cf*Lf-Cr*Lr)/(m*Vx); 
    0 0 0 1; 
    0 -2*(Cf*Lf-Cr*Lr)/(Iz*Vx) 2*(Cf*Lf-Cr*Lr)/Iz -2*(Cf*Lf^2+Cr*Lr^2)/(Iz*Vx)];
B=[0; -Vx-2*(Cf*Lf-Cr*Lr)/(m*Vx); 0 ;-2*(Cf*Lf^2+Cr*Lr^2)/(Iz*Vx)];
Bv = [0;2*Cf/m;0;2*Cf*Lf/Iz];
D = zeros(4,2);

Q=[.8 0 0 0; 0 .5 0 0; 0 0 .2 0; 0 0 0 .6];
R=0.4;
[X,K,L] = idare(A,B,Q,R,[],[])
C=[1 0 0 -1; 0 1 0 0];
ns=size(A,1);
na=size(B,2); 
nm=size(C,1); 
nd=na+nm;     
nr=nm+na;     

B1=[B zeros(4)];
B2=[B];
C1=[C; zeros(4,4)];
C2=[C];
D11=zeros(6,4);
D12=[D; eye(2)];
D21=[D, eye(4,2)];
D22=D;
%%
Fd = sdpvar(1,4);
P=sdpvar(4,4);
gamma=sdpvar(1);
F=[];
F1=[P A*P-B2*F B1 zeros(4,1);
    (A*P-B2*F)' P  [] (P*C1'-F'*D11'); 
    B1' zeros(1,4) gamma 0;
    zeros(4,2) (P*C1'-F*D11')' 0 gamma];
F=[F;F1>0];
F=[F;P>0.001*eye(4)];

optimize(F,gamma)

value(P)

%%
%Formulating the state space of the 9-matrix rep to find Hinf-norm
sys=ss(A ,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
norm(sys,inf)

%%
F=[];
F1=[A'*P*A-P+Q A'*P*B; 
    B'*P*A R+B'*P*B]
F=[F;F1>=0];
F=[F;P>0.001*eye(4)];

optimize(F,[])

value(P)

%% DARE
P=sdpvar(4,4);
Q=[.8 0 0 0; 0 .5 0 0; 0 0 .2 0; 0 0 0 .6];
R=0.4;

J = [A'*P*A-P+Q A'*P*B;(A'*P*B)' R+B'*P*B];

F=[P>1e-6*eye(4)];
F=[F;J>=0];

optimize(F,P)
value(P)

 