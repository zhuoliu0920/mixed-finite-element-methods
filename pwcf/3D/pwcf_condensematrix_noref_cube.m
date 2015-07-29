
h = 1;

M1 = 1/6*h*h.*[3,sqrt(3)/2,sqrt(3)/2,sqrt(3);sqrt(3)/2,1,0,1/2;sqrt(3)/2,0,1,1/2;sqrt(3),1/2,1/2,2];
B1 = -h*h/2.*[sqrt(3),1,1,1];
C1 = diag(-B1);
invM1 = inv(M1);
b1 = B1*inv(M1)*B1';
A1 = inv(M1)*B1';
D1 = (1/b1).*(A1*A1');
E1 = inv(M1)-D1;
S1 = C1*E1*C1';