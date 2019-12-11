A= [-0.0558 -0.9968 0.0802 0.0415
0.5980 -0.1150 -0.0318 0
-3.0500 0.3880 -0.4650 0
0 0.0805 1.0000 0];

B= [0.0729 0.0001
-4.7500 1.2300
1.5300 10.6300
0 0];

C= [ 0 1 0 0; 0 0 0 1];

P = sdpvar(4,4);
Z = sdpvar(2,2);
V= sdpvar(4,4);
gamma = sdpvar(1);

mat1 = [-(V+V') V'*A'+P V'*C' V';
A*V+P -P zeros(4,2) zeros(4,4);
C*V zeros(2,4) -eye(2) zeros(2,4);
V zeros(4) zeros(4,2) -P];

mat2 = [-Z B'; B -P];

F = [mat1 < 0; mat2<0; P>0; trace(Z)<=gamma];

optimize(F, gamma);

H2_norm = sqrt(value(gamma))