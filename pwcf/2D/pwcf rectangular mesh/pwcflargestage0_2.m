function [u,p,lambda] = pwcflargestage0_2(n)
stage = 1;
h = 1/n;
dimM = 4*(n*n+3*stage);
dimB = n*n+3*stage;
dimC = 2*n*(n-1)+6*stage;

% generating matrices Mi and Bi
subM = sparse((h*h/2)*eye(4));
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
l1 = zeros(n*(n-2),1);
m1 = zeros(n*(n-2),1);
l2 = zeros(n*(n-2),1);
m2 = zeros(n*(n-2),1);
for j = 1:(n-2)
    for i = 1:n
        l1(i+(j-1)*n) = i+(j-1)*n;
        m1(i+(j-1)*n) = (i-1)*4*n+1+(j-1)*4;
        l2(i+(j-1)*n) = i+(j-1)*n;
        m2(i+(j-1)*n) = (i-1)*4*n+7+(j-1)*4;
    end
end
    j = n-1;
    for i = 1:n-1
        l1(i+(j-1)*n) = i+(j-1)*n;
        m1(i+(j-1)*n) = (i-1)*4*n+1+(j-1)*4;
        l2(i+(j-1)*n) = i+(j-1)*n;
        m2(i+(j-1)*n) = (i-1)*4*n+7+(j-1)*4;
    end
s1 = h*ones(n*(n-1)-1,1);
s2 = -s1;
C = sparse(l1,m1,s1,dimC,dimM) + sparse(l2,m2,s2,dimC,dimM);

% the vertical intersections
for k = 1:(n*(n-1)-1)
    l1(k) = k+n*(n-1)-1;
    m1(k) = (k-1)*4+2;
    l2(k) = k+n*(n-1)-1;
    m2(k) = (k-1)*4+4+4*n;
end
C = C + sparse(l1,m1,s1,dimC,dimM) + sparse(l2,m2,s2,dimC,dimM);
p = 2*(n*(n-1)-1);

C(p+1,4*(n*n-1)-3) = h/2;
C(p+1,4*n*n-1) = -h/2;
C(p+2,4*(n*n-1)-3) = h/2;
C(p+2,4*(n*n+2)-1) = -h/2;
C(p+3,4*n*(n-1)-2) = h/2;
C(p+3,4*n*n) = -h/2;
C(p+4,4*n*(n-1)-2) = h/2;
C(p+4,4*(n*n+1)) = -h/2;
C(p+5,4*n*n-3) = h/2;
C(p+5,4*(n*n+1)-1) = -h/2;
C(p+6,4*(n*n+2)-3) = h/2;
C(p+6,4*(n*n+3)-1) = -h/2;
C(p+7,4*(n*n)-2) = h/2;
C(p+7,4*(n*n+2)) = -h/2;
C(p+8,4*(n*n+1)-2) = h/2;
C(p+8,4*(n*n+3)) = -h/2;

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
for s = 1:stage
    F(n*n+3*(s-1)) = -2*(cos(pi*(1/2)^(s-1)/n)-cos(pi*(1/2)^(s)/n))^2;
    F(n*n+3*(s-1)+1) = -2*(cos(pi*(1/2)^(s-1)/n)-cos(pi*(1/2)^(s)/n))*(cos(pi*(1/2)^(s)/n)-1);
    F(n*n+3*(s-1)+2) = -2*(cos(pi*(1/2)^(s)/n)-1)*(cos(pi*(1/2)^(s-1)/n)-cos(pi*(1/2)^(s)/n));
    F(n*n+3*(s-1)+3) = -2*(cos(pi*(1/2)^(s)/n)-1)^2;
end

% generating the matrix S and vector psi
S = sparse(zeros(dimC,dimC));
psi = zeros(dimC,1);
for i = 1:(n*n-1)
    subC = C(:,(4*(i-1)+1):(4*i));
    S = S + subC*inv(subM)*(subM-subB'*inv(subB*inv(subM)*subB')*subB)*inv(subM)*subC';
    psi = psi + subC*(inv(subM)*subB'*inv(subB*inv(subM)*subB'))*F(i);
end

for s = 1:stage
    subM = subM/4;
    subB = subB/2;
    for j = 0:2
        i = n*n+3*(s-1)+j;
        subC = C(:,(4*(i-1)+1):(4*i));
        S = S + subC*inv(subM)*(subM-subB'*inv(subB*inv(subM)*subB')*subB)*inv(subM)*subC';
        psi = psi + subC*(inv(subM)*subB'*inv(subB*inv(subM)*subB'))*F(i);
    end
end
subC = C(:,(dimM-3):dimM);
S = S + subC*inv(subM)*(subM-subB'*inv(subB*inv(subM)*subB')*subB)*inv(subM)*subC';
psi = psi + subC*(inv(subM)*subB'*inv(subB*inv(subM)*subB'))*F(n*n+3*stage);
lambda = S\psi;

 
%redefine matrices subM and subB: 
subM = sparse((h*h/2)*eye(4));
subB = zeros(1,4);
  for  i = 1:2
     subB(i) = -h;
  end   
  for  i = 3:4
     subB(i) = h; 
  end
subB = sparse(subB);
% computing p
p = zeros(dimB,1);
for i = 1:(n*n-1)
    subC = C(:,(4*(i-1)+1):(4*i));
    S2 = subB*inv(subM)*subB';
    psi2 = -subB*inv(subM)*subC'*lambda-F(i);
    p(i) = psi2/S2;
end
for s = 1:stage
    subM = subM/4;
    subB = subB/2;
    for j = 0:2
        i = n*n+3*(s-1)+j;
        subC = C(:,(4*(i-1)+1):(4*i));
        S2 = subB*inv(subM)*subB';
        psi2 = -subB*inv(subM)*subC'*lambda-F(i);
        p(i) = psi2/S2;
    end
end
subC = C(:,(dimM-3):dimM);
S2 = subB*inv(subM)*subB';
psi2 = -subB*inv(subM)*subC'*lambda-F(n*n+3*stage);
p(n*n+3*stage) = psi2/S2;

%redefine matrices subM and subB: 
subM = sparse((h*h/2)*eye(4));
subB = zeros(1,4);
  for  i = 1:2
     subB(i) = -h;
  end   
  for  i = 3:4
     subB(i) = h; 
  end
subB = sparse(subB);
% computing u
u = zeros(dimM,1);
for i = 1:(n*n-1)
    subC = C(:,(4*(i-1)+1):(4*i));
    psi3 = -subC'*lambda-subB'*p(i);
    u(4*(i-1)+1:4*i) = subM\psi3;
end
for s = 1:stage
    subM = subM/4;
    subB = subB/2;
    for j = 0:2
        i = n*n+3*(s-1)+j;
        subC = C(:,(4*(i-1)+1):(4*i));
        psi3 = -subC'*lambda-subB'*p(i);
        u(4*(i-1)+1:4*i) = subM\psi3;
    end
end
subC = C(:,(dimM-3):dimM);
S2 = subB*inv(subM)*subB';
psi3 = -subC'*lambda-subB'*p(n*n+3*stage);
u(4*(n*n+3*stage-1)+1:4*(n*n+3*stage)) = subM\psi3;

