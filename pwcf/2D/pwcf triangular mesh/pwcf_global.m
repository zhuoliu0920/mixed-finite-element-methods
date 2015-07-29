A = sqrt(3);
a = A/4;
Msub = [2/3,1/3,0;1/3,4/3,1/3;0,1/3,2/3];
M1 = A*Msub;
M2 = [8*a/3,a/3,a/3,a/3,a/3;a/3,2*a/3,0,0,0;a/3,0,2*a/3,0,0;a/3,0,0,2*a/3,0;a/3,0,0,0,2*a/3];
M3 = a*Msub;
M4 = a*Msub;
M5 = a*Msub;
M = blkdiag(M1,M2,M3,M4,M5);

S = 2;
s = 1;
B1 = [S,S,S];
B2 = [s,s,s,0,0;s,0,0,s,s];
B3 = [s,s,s];
B4 = [s,s,s];
B5 = [s,s,s];
B = -blkdiag(B1,B2,B3,B4,B5);

Ctrans = [0,0,0,0,0;S,0,0,0,0;0,0,0,0,0;S,0,0,0,0;0,s,0,0,0;0,0,s,0,0;0,0,0,s,0;0,0,0,0,s;0,0,0,0,0;0,0,0,0,0;0,s,0,0,0;0,0,s,0,0;0,0,0,0,0;0,0,0,s,0;0,0,0,0,s;0,0,0,0,0;0,0,0,0,0];
C = Ctrans';

Global = [M B' C';B zeros(6,11);C zeros(5,11)];
F = [zeros(17,1);18;9/2;9/2;9/2;9/2;9/2;zeros(5,1)];
f = [18;9/2;9/2;9/2;9/2;9/2];
sol = linsolve(Global,F);

J = B*inv(M)*C'*inv(C*inv(M)*C')*C*inv(M)*B';
H = -(B*inv(M)*B'-B*inv(M)*C'*inv(C*inv(M)*C')*C*inv(M)*B');
solp = linsolve(H,f);

%[X,Y] = meshgrid(0:.02:3, 0:.02:3*sqrt(3)/2);                                
%Z = Y.*(Y-sqrt(3).*X).*(3.*sqrt(3)-sqrt(3).*X-Y);                                     
%surf(X,Y,Z);