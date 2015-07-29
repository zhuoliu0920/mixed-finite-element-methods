function [u,p,lambda, M, B, C, H,H1,H2, Global] = rt0_tri_noref(n)

h1 = sqrt(3)/n;
b1 = 2/n;
a1 = h1*b1/2;

dimM = 3*n*n;
dimB = n*n;
dimC = 3*n*(n-1)/2;

Mtri = (1/3).*[5/3,-1/3,-1/3;-1/3,5/3,-1/3;-1/3,-1/3,5/3];
M = [];
for i = 1:n*n
    M = blkdiag(M,a1.*Mtri);
end

Btri1 = -[1,1,1];
Btri2 = [1,1,1];
B = [];

row_ele(1) = 1;
total_ele(1) = 1;
for i = 2:n
    row_ele(i) = 2*i-1;
    total_ele(i) = total_ele(i-1)+row_ele(i);
end

for i = 1:n
    for j = 1:row_ele(i)
        if rem(j,2) == 0
            B = blkdiag(B, b1.*Btri2);
        else
            B = blkdiag(B, b1.*Btri1);
        end
    end
end

C = zeros(dimC,dimM);
seq1(1) = 3;
row1(1) = 1;
for i = 2:(n-1)
    seq1((i*i-i)/2+1) = seq1((i*i-i)/2)+3;
    row1((i*i-i)/2+1) = row1((i*i-i)/2)+1;
    for j = 2:i
        seq1((i*i-i)/2+j) = seq1((i*i-i)/2+j-1)+2;
        row1((i*i-i)/2+j) = row1((i*i-i)/2+j-1);
    end
end

for i = 1:(n*n/2-n/2)
    C(3*i-2,3*seq1(i)-1) = -b1; C(3*i-2,3*(seq1(i)-2*row1(i))-1) = b1; 
    C(3*i-1,3*seq1(i)-2) = -b1; C(3*i-1,3*(seq1(i)-1)) = b1;
    C(3*i,3*seq1(i)) = -b1; C(3*i,3*seq1(i)+1) = b1;
end

f = zeros(n*n, 1);
for i = 1:n*n
    f(i) = a1*(-4*sqrt(3));
end

F = [zeros(dimM,1);f;zeros(dimC,1)];
Global = [M B' C';B zeros(dimB,dimB+dimC);C zeros(dimC,dimB+dimC)];
sol = linsolve(Global,F);
u = sol(1:dimM); p = sol((dimM+1):(dimM+dimB)); lambda = sol((dimM+dimB+1):(dimM+dimB+dimC));

%F = [zeros(17,1);18;9/2;9/2;9/2;9/2;9/2;zeros(5,1)];
%f = [18;9/2;9/2;9/2;9/2;9/2];
%sol = linsolve(Global,F);

%J = B*inv(M)*C'*inv(C*inv(M)*C')*C*inv(M)*B';
H = (B*inv(M)*B'-B*inv(M)*C'*inv(C*inv(M)*C')*C*inv(M)*B');
H1 = C*inv(M)*(M-B'*inv(B*inv(M)*B')*B)*inv(M)*C';
h11 = B*inv(M)*B'; h12 = B*inv(M)*C'; h21 = C*inv(M)*B'; h22 = C*inv(M)*C';
H2 = -[h11 h12; h21 h22];
%solp = linsolve(H,f);

%[X,Y] = meshgrid(0:.02:3, 0:.02:3*sqrt(3)/2);                                
%Z = Y.*(Y-sqrt(3).*X).*(3.*sqrt(3)-sqrt(3).*X-Y);                                     
%surf(X,Y,Z);
