function [u,p,lambda, M, B, C, H, S, Global] = pwcf_tri_ref(n)

h1 = sqrt(3)/(2*n);
h2 = h1/2;
b1 = 1/n;
b2 = b1/2;
a1 = h1*b1/2;
a2 = h2*b2/2;

dimM = 39*n*n-n;
dimB = 13*n*n;
dimC = (3*n*n-3*n)/2+5*n+18*n*n-9*n;

Mtri = [2/3,1/3,0;1/3,4/3,1/3;0,1/3,2/3];
Mpen = [8/3,1/3,1/3,1/3,1/3;1/3,2/3,0,0,0;1/3,0,2/3,0,0;1/3,0,0,2/3,0;1/3,0,0,0,2/3];
M = [];
for i = 1:n*n
    M = blkdiag(M,a1.*Mtri);
end
for i = (n*n+1):(n*n+n)
    M = blkdiag(M,a2.*Mpen);
end
for i = (n*n+n+1):(13*n*n-n)
    M = blkdiag(M,a2.*Mtri);
end

Btri1 = -[1,1,1];
Btri2 = [1,1,1];
Bpen = [1,1,1,0,0;1,0,0,1,1];
B = [];

row_ele(1) = 1;
total_ele(1) = 1;
for i = 2:n
    row_ele(i) = 2*i-1;
    total_ele(i) = total_ele(i-1)+row_ele(i);
end
    total_ele(n+1) = total_ele(n)+3*n+1;
for i = (n+2):3*n
    row_ele(i) = 4*n-1+2*(i-n);
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
for i = 1:n
    B = blkdiag(B, b2.*Bpen);
end
for i = 1:(2*n+1)
    B = blkdiag(B, b2.*Btri1);
end
for i = (n+2):3*n
    for j = 1:row_ele(i)
        if rem(j,2) == 0
            B = blkdiag(B, b2.*Btri2);
        else
            B = blkdiag(B, b2.*Btri1);
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
for i = 1:n
    C(3*n*n/2-3*n/2+5*i-4,3*(n*n-2*n+2*i)-1) = b1; C(3*n*n/2-3*n/2+5*i-4,3*n*n+5*(i-1)+1) = -b1;
    C(3*n*n/2-3*n/2+5*i-3,3*(n*n+2*i-1)+5*n) = b2; C(3*n*n/2-3*n/2+5*i-3,3*n*n+5*(i-1)+2) = -b2;
    C(3*n*n/2-3*n/2+5*i-2,3*(n*n+2*i-1)+5*n+1) = b2; C(3*n*n/2-3*n/2+5*i-2,3*n*n+5*(i-1)+3) = -b2;
    C(3*n*n/2-3*n/2+5*i-1,3*(n*n+2*i-1)+5*n+3) = b2; C(3*n*n/2-3*n/2+5*i-1,3*n*n+5*(i-1)+4) = -b2;
    C(3*n*n/2-3*n/2+5*i,3*(n*n+2*i-1)+5*n+4) = b2;C(3*n*n/2-3*n/2+5*i,3*n*n+5*(i-1)+5) = -b2;
end

seq2(1) = 2;
row2(1) = 1;
pos2(1) = 2*n+2;
for i = 1:(2*n-1)
    for j = 2:(2*n+i)
        seq2(2*n*(i-1)+i*(i-1)/2+j) = seq2(2*n*(i-1)+i*(i-1)/2+j-1)+2;
        row2(2*n*(i-1)+i*(i-1)/2+j) = row2(2*n*(i-1)+i*(i-1)/2+j-1);
        if i == 1
            pos2(2*n*(i-1)+i*(i-1)/2+j) = pos2(2*n*(i-1)+i*(i-1)/2+j-1)+1;
        else
            pos2(2*n*(i-1)+i*(i-1)/2+j) = pos2(2*n*(i-1)+i*(i-1)/2+j-1);
        end
    end
    seq2(2*n*i+i*(i+1)/2+1) = seq2(2*n*i+i*(i+1)/2)+3;
    row2(2*n*i+i*(i+1)/2+1) = row2(2*n*i+i*(i+1)/2)+1;
    pos2(2*n*i+i*(i+1)/2+1) = pos2(2*n*i+i*(i+1)/2)+2;
end
for i = 1:(6*n*n-3*n)
    C(3*n*n/2-3*n/2+5*n+3*i-2,3*n*n+11*n+3+3*seq2(i)-1) = -b2; C(3*n*n/2-3*n/2+5*n+3*i-2,3*n*n+11*n+3+3*(seq2(i)-pos2(i))-1) = b2;
    C(3*n*n/2-3*n/2+5*n+3*i-1,3*n*n+11*n+3+3*seq2(i)-2) = -b2; C(3*n*n/2-3*n/2+5*n+3*i-1,3*n*n+11*n+3+3*seq2(i)-3) = b2;
    C(3*n*n/2-3*n/2+5*n+3*i,3*n*n+11*n+3+3*seq2(i)) = -b2; C(3*n*n/2-3*n/2+5*n+3*i,3*n*n+11*n+3+3*seq2(i)+1) = b2;
end

f = zeros(13*n*n, 1);
for i = 1:n*n
    f(i) = a1*(-4*sqrt(3));
end
for i = (n*n+1):(13*n*n)
    f(i) = a2*(-4*sqrt(3));
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
%solp = linsolve(H,f);
S = C*(inv(M)-inv(M)*B'*inv(B*inv(M)*B')*B*inv(M))*C';

%[X,Y] = meshgrid(0:.02:3, 0:.02:3*sqrt(3)/2);                                
%Z = Y.*(Y-sqrt(3).*X).*(3.*sqrt(3)-sqrt(3).*X-Y);                                     
%surf(X,Y,Z);
