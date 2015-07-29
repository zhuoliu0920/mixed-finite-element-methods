clear;
n = 1;
h1 = sqrt(3)/(2*n);
h2 = h1/2;
b1 = 1/(2*n-1);
b2 = b1/2;
a1 = h1*b1/2;
a2 = h2*b2/2;

dimM = 39*n*n-n;
dimB = 13*n*n;
dimC = (3*n*n-3*n)/2+5*n+18*n*n-9*n;

Mtri = [2/3,1/3,0;1/3,4/3,1/3;0,1/3,2/3];
Mpen = [8/3,1/3,1/3,1/3,1/3;1/3,2/3,0,0,0;1/3,0,2/3,0,0;1/3,0,0,2/3,0;1/3,0,0,0,2/3];
Btri = -[1,1,1];
Bpen = -[1,1,1,0,0;1,0,0,1,1];

M = [];
B = [];
for i = 1:n*n
    M = blkdiag(M,a1.*Mtri);
    B = blkdiag(B,b1.*Btri);
end
for i = (n*n+1):(n*n+n)
    M = blkdiag(M,a2.*Mpen);
    B = blkdiag(B,b2.*Bpen);
end
for i = (n*n+n+1):(13*n*n-n)
    M = blkdiag(M,a2.*Mtri);
    B = blkdiag(B,b2.*Btri);
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
    C(3*i-2,3*seq1(i)-1) = b1; C(3*i-2,3*(seq1(i)-2*row1(i))-1) = b1; 
    C(3*i-1,3*seq1(i)-2) = b1; C(3*i-1,3*(seq1(i)-1)) = b1;
    C(3*i,3*seq1(i)) = b1; C(3*i,3*seq1(i)+1) = b1;
end
for i = 1:n
    C(3*n*n/2-3*n/2+5*i-4,3*(n*n-2*n+2*i)-1) = b1; C(3*n*n/2-3*n/2+5*i-4,3*(n*n+i-1)+1) = b1;
    C(3*n*n/2-3*n/2+5*i-3,3*(n*n+2*i-1)+5*n) = b2; C(3*n*n/2-3*n/2+5*i-3,3*(n*n+i-1)+2) = b2;
    C(3*n*n/2-3*n/2+5*i-2,3*(n*n+2*i-1)+5*n+1) = b2; C(3*n*n/2-3*n/2+5*i-2,3*(n*n+i-1)+3) = b2;
    C(3*n*n/2-3*n/2+5*i-1,3*(n*n+2*i-1)+5*n+3) = b2; C(3*n*n/2-3*n/2+5*i-1,3*(n*n+i-1)+4) = b2;
    C(3*n*n/2-3*n/2+5*i,3*(n*n+2*i-1)+5*n+4) = b2;C(3*n*n/2-3*n/2+5*i,3*(n*n+i-1)+5) = b2;
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
    C(3*n*n/2-3*n/2+5*n+3*i-2,3*n*n+11*n+3+3*seq2(i)-1) = b2; C(3*n*n/2-3*n/2+5*n+3*i-2,3*n*n+11*n+3+3*(seq2(i)-pos2(i))-1) = b2;
    C(3*n*n/2-3*n/2+5*n+3*i-1,3*n*n+11*n+3+3*seq2(i)-2) = b2; C(3*n*n/2-3*n/2+5*n+3*i-1,3*n*n+11*n+3+3*seq2(i)-3) = b2;
    C(3*n*n/2-3*n/2+5*n+3*i,3*n*n+11*n+3+3*seq2(i)) = b2; C(3*n*n/2-3*n/2+5*n+3*i,3*n*n+11*n+3+3*seq2(i)+1) = b2;
end

%Global = [M B' C';B zeros(6,11);C zeros(5,11)];
%F = [zeros(17,1);18;9/2;9/2;9/2;9/2;9/2;zeros(5,1)];
%f = [18;9/2;9/2;9/2;9/2;9/2];
%sol = linsolve(Global,F);

%J = B*inv(M)*C'*inv(C*inv(M)*C')*C*inv(M)*B';
H = (B*inv(M)*B'-B*inv(M)*C'*inv(C*inv(M)*C')*C*inv(M)*B');
%solp = linsolve(H,f);

%[X,Y] = meshgrid(0:.02:3, 0:.02:3*sqrt(3)/2);                                
%Z = Y.*(Y-sqrt(3).*X).*(3.*sqrt(3)-sqrt(3).*X-Y);                                     
%surf(X,Y,Z);
