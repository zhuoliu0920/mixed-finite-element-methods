clear;
% input coordinates of four vertexes.

h = 1;

%x1 = [0,0,1]; x4 = [0,0,0]; x2 = [0,1,0]; x3 = [1,0,0];

%x1 = h.*[sqrt(3)/3,0,2*sqrt(6)/3]; x2 = h.*[0,-1,0]; x3 = h.*[0,1,0]; x4 = h.*[sqrt(3),0,0];

x1 = [0,0,0.01]; x2 = [-1,0,0]; x3 = [1,0,0]; x4 = [0,1,0]; 

vol = 1/6*abs(det([x2-x1;x3-x1;x4-x1]));

% step 1: generating normal fluxes on each interfaces.
n1 = cross(x2-x4,x3-x4); s1 = norm(n1)/2; n1 = n1/norm(n1); 
n2 = cross(x3-x4,x1-x4); s2 = norm(n2)/2; n2 = n2/norm(n2);
n3 = cross(x2-x1,x4-x1); s3 = norm(n3)/2; n3 = n3/norm(n3);
n4 = cross(x3-x1,x2-x1); s4 = norm(n4)/2; n4 = n4/norm(n4);

n12 = dot(n1,n2); n13 = dot(n1,n3); n14 = dot(n1,n4); n23 = dot(n2,n3); n24 = dot(n2,n4); n34 = dot(n3,n4);

% step 2: finding basis for PWC fluxes on each interfaces. 
w11 = linsolve([n1;n3;n4],[1;0;0]); w12 = linsolve([n1;n2;n4],[1;0;0]);
w2 = linsolve([n2;n1;n4],[1;0;0]);
w3 = linsolve([n3;n1;n4],[1;0;0]);
w41 = linsolve([n4;n3;n1],[1;0;0]); w42 = linsolve([n4;n2;n1],[1;0;0]);

% step 3: generating local MHFE matrices.

M = vol/2.*[dot(w11,w11)+dot(w12,w12),dot(w12,w2),dot(w11,w3),dot(w11,w41)+dot(w12,w42);
            dot(w2,w12),dot(w2,w2),0,dot(w2,w42);
            dot(w3,w11),0,dot(w3,w3),dot(w3,w41);
            dot(w41,w11)+dot(w42,w12),dot(w42,w2),dot(w41,w3),dot(w41,w41)+dot(w42,w42)];
condPWCF = cond(M);