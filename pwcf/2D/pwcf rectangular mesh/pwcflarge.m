function [u,p,lambda] = pwcflarge(n)
h = 1/n;
dimM = 4*n*n;
dimB = n*n;
dimC = 2*n*(n-1);

% generating matrices Mi and Bi
subM = sparse(h*h/2*eye(4));
subB = zeros(1,4);
  for  i = 1:2
     subB(i) = -h;
  end   
  for  i = 3:4
     subB(i) = h; 
  end
subB = sparse(subB);

% generating the sparse matrix C:
% the horizontal intersections
l1 = zeros(n*(n-1),1);
m1 = zeros(n*(n-1),1);
l2 = zeros(n*(n-1),1);
m2 = zeros(n*(n-1),1);
for i = 1:n
    for j = 1:(n-1)
        l1(i+(j-1)*n) = i+(j-1)*n;
        m1(i+(j-1)*n) = (i-1)*4*n+1+(j-1)*4;
        l2(i+(j-1)*n) = i+(j-1)*n;
        m2(i+(j-1)*n) = (i-1)*4*n+7+(j-1)*4;
    end
end
s1 = h*ones(n*(n-1),1);
s2 = -s1;
C = sparse(l1,m1,s1,dimC,dimM)+sparse(l2,m2,s2,dimC,dimM);
% the vertical intersections
for k = 1:(n*(n-1))
    l1(k) = k+n*(n-1);
    m1(k) = (k-1)*4+2;
    l2(k) = k+n*(n-1);
    m2(k) = (k-1)*4+4+4*n;
end
C = C + sparse(l1,m1,s1,dimC,dimM)+sparse(l2,m2,s2,dimC,dimM);

% generating F
% delta p = f in omega;
%       p = 0 on boundary.
% so p = sin(pi*x1)*sin(pi*x2), f = -2*pi*pi*sin(pi*x1)*sin(pi*x2)
F = zeros(dimB,1);
for i = 0:(n-1)
    for j = 0:(n-1)
    F(n*i+j+1) = -2*(cos(pi-i*pi/n)-cos(pi-(i+1)*pi/n))*(cos(pi-j*pi/n)-cos(pi-(j+1)*pi/n));
    end
end
% generating the matrix S and vector psi
S = sparse(zeros(dimC,dimC));
psi = zeros(dimC,1);
for i = 1:dimB
    subC = C(:,(4*(i-1)+1):(4*i));
    S = S + subC*inv(subM)*(subM-subB'*inv(subB*inv(subM)*subB')*subB)*inv(subM)*subC';
    psi = psi + subC*(inv(subM)*subB'*inv(subB*inv(subM)*subB'))*F(i);
end
lambda = S\psi;

% computing p
p = zeros(dimB,1);
for i = 1:dimB
    subC = C(:,(4*(i-1)+1):(4*i));
    S2 = subB*inv(subM)*subB';
    psi2 = -subB*inv(subM)*subC'*lambda-F(i);
    p(i) = psi2/S2;
end

% computing u
u = zeros(dimM,1);
for i = 1:dimB
    subC = C(:,(4*(i-1)+1):(4*i));
    psi3 = -subC'*lambda-subB'*p(i);
    u(4*(i-1)+1:4*i) = subM\psi3;
end
