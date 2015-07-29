a1 = pi/3; a2 = pi/3;
k = 1/2;
M1 = [1/sin(a1)^2,cos(a1)/sin(a1)^2,0;cos(a1)/sin(a1)^2,1/sin(a1)^2,0;0,0,0];
M2 = [0,0,0;0,1/sin(a2)^2,cos(a2)/sin(a2)^2;0,cos(a2)/sin(a2)^2,1/sin(a2)^2];
M = (cot(a1)/2+cot(a2)/2).*(k.*M1+(1-k).*M2);
B = [-1/sin(a1),-cot(a1)-cot(a2),-1/sin(a2)];
C = diag(-B);
invM = inv(M);
b = B*inv(M)*B';
A = inv(M)*B';
D = (1/b).*(A*A');
E = inv(M)-D;
S = C*E*C';
cot1 = cot(a1);
cot2 = cot(a2);

S1 = (1/(cot1/2+cot2/2)).*[csc(a1)^2,-cot1^2-cot1*cot2,cot1*cot2-1;-cot1^2-cot1*cot2,(cot1+cot2)^2,-cot2^2-cot1*cot2;cot1*cot2-1,-cot2^2-cot1*cot2,csc(a2)^2];