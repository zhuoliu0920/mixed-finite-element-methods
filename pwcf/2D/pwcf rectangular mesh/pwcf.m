n = 3;
h = 1/n;
dimM = 4*n*n;
dimB = n*n;
dimC = 2*n*(n-1);
M = h*h/2*eye(dimM);
B = zeros(dimB,dimM);
for i = 1:dimB
  for  j = 1:2
     B(i,(i-1)*4+j) = -h;
  end   
  for  j = 3:4
     B(i,(i-1)*4+j) = h; 
  end
end
C = zeros(dimC,dimM);
% the horizontal intersections
for i = 1:n
    for j = 1:(n-1)
        C(i+(j-1)*n, (i-1)*4*n+1+(j-1)*4) = 1;
        C(i+(j-1)*n, (i-1)*4*n+7+(j-1)*4) = -1;
    end
end
% the vertical intersections
for k = 1:(n*(n-1))
        C(k+n*(n-1),(k-1)*4+2) = 1;
        C(k+n*(n-1),(k-1)*4+4+4*n) = -1;
end
C = h*C;
clear i j k;

% delta p = f in omega;
%       p = 0 on boundary.
% so p = sin(pi*x1)*sin(pi*x2), f = -2*pi*pi*sin(pi*x1)*sin(pi*x2)
F = zeros(dimB,1);
for i = 0:(n-1)
    for j = 0:(n-1)
    F(n*i+j+1) = -2*(cos(pi-i*pi/n)-cos(pi-(i+1)*pi/n))*(cos(pi-j*pi/n)-cos(pi-(j+1)*pi/n));
    end
end
clear i j;
A = [M, B', C'; B, zeros(dimB,dimB), zeros(dimB,dimC); C, zeros(dimC,dimB), zeros(dimC,dimC)];
f = [zeros(dimM,1); F; zeros(dimC,1)];
U = linsolve(A,f);
u = U(1:dimM);
p = U((dimM+1):(dimM+dimB));
lambda = U((dimM+dimB+1):(dimM+dimB+dimC));
clear dimM dimB dimC
