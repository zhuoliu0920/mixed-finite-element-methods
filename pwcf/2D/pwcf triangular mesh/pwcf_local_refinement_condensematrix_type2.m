clear;
a1 = pi/3; a2 = pi/3;

M = [cos(a1)/sin(a1)^3,cos(a1)^2/sin(a1)^3,0,0,0; cos(a1)^2/sin(a1)^3,cos(a1)/sin(a1)^3+cos(a2)/sin(a2)^3+cot(a1)+cot(a2),0,0,cos(a2)^2/sin(a2)^3; 0,0,cot(a1),0,0;0,0,0,cot(a2),0; 0,cos(a2)^2/sin(a2)^3,0,0,cos(a2)/sin(a2)^3];

B = [-1/sin(a1),-cos(a1)/sin(a1),-1,0,0; 0,-cos(a2)/sin(a2),0,1,-1/sin(a2)];
C = zeros(4,5);
C(1,1) = 1/sin(a1); C(2,2) = cot(a1)+cot(a2); C(3,3) = 1; C(3,4) = -1; C(4,5) = 1/sin(a2);

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