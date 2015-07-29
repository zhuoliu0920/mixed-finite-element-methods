clear;
a1 = pi/3; a2 = pi/3; a3 = pi/3; a4 = pi/3;
k1 = 1/2; k2 = 1/2;
M1 = [1/sin(a1)^2,cos(a1)/sin(a1)^2,0,0,0;cos(a1)/sin(a1)^2,1/sin(a1)^2,0,0,0;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0];
M2 = [1/sin(a2)^2,0,cos(a2)/sin(a2)^2,0,0;0,0,0,0,0;cos(a2)/sin(a2)^2,0,1/sin(a2)^2,0,0;0,0,0,0,0;0,0,0,0,0];
M3 = [1/sin(a3)^2,0,0,cos(a3)/sin(a3)^2,0;0,0,0,0,0;0,0,0,0,0;cos(a3)/sin(a3)^2,0,0,1/sin(a3)^2,0;0,0,0,0,0];
M4 = [1/sin(a4)^2,0,0,0,cos(a4)/sin(a4)^2;0,0,0,0,0;0,0,0,0,0;0,0,0,0,0;cos(a4)/sin(a4)^2,0,0,0,1/sin(a4)^2];

M = (k1.*M1+(1-k1).*M2)+(k2.*M3+(1-k2).*M4);
%M = k1.*M1+(1-k1).*M2+k2.*M3+(1-k2).*M4;
B = [-cot(a1)-cot(a2),-1/sin(a1),-1/sin(a2),0,0;-cot(a3)-cot(a4),0,0,-1/sin(a3),-1/sin(a4)];
B1 = [-cot(a1)-cot(a2)-cot(a3)-cot(a4),-1/sin(a1),-1/sin(a2),-1/sin(a3),-1/sin(a4)];
C = diag(-B1);
invM = inv(M);
b = B*inv(M)*B';
A = inv(M)*B';
D = A*inv(b)*A';
E = inv(M)-D;
S = C*E*C';
% cot1 = cot(a1);
% cot2 = cot(a2);

% invM2 = det(M)*(sin(a1)*sin(a2)*sin(a3)*sin(a4))^2/(k1*(1-k1)*k2*(1-k2)).*invM;

% S1 = (1/(cot1/2+cot2/2)).*[csc(a1)^2,-cot1^2-cot1*cot2,cot1*cot2-1;-cot1^2-cot1*cot2,(cot1+cot2)^2,-cot2^2-cot1*cot2;cot1*cot2-1,-cot2^2-cot1*cot2,csc(a2)^2];