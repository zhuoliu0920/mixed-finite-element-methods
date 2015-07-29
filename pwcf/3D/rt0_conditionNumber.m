clear;
% input coordinates of four vertexes.

h = 1;

%x1 = [0,0,1]; x4 = [0,0,0]; x2 = [0,1,0]; x3 = [1,0,0];
 
%x1 = h.*[sqrt(3)/3,0,2*sqrt(6)/3]; x2 = h.*[0,-1,0]; x3 = h.*[0,1,0]; x4 = h.*[sqrt(3),0,0];

x1 = [0,0,0.01]; x2 = [-1,0,0]; x3 = [1,0,0]; x4 = [0,1,0]; 

vol = 1/6*abs(det([x2-x1;x3-x1;x4-x1]));

% step 1: generating normal fluxes on each interfaces.
n1 = cross(x2-x4,x3-x4); s1 = norm(n1)/2; n1 = n1/norm(n1); h1 = ptPlaneDist(x1, x2, x3, x4); 
n2 = cross(x3-x4,x1-x4); s2 = norm(n2)/2; n2 = n2/norm(n2); h2 = ptPlaneDist(x2, x1, x3, x4);
n3 = cross(x2-x1,x4-x1); s3 = norm(n3)/2; n3 = n3/norm(n3); h3 = ptPlaneDist(x3, x1, x2, x4);
n4 = cross(x3-x1,x2-x1); s4 = norm(n4)/2; n4 = n4/norm(n4); h4 = ptPlaneDist(x4, x1, x2, x3);

% step 2: finding RT0 basis.
syms x y z;
P = [x,y,z];
w1 = (P-x1)./h1; 
w2 = (P-x2)./h2; 
w3 = (P-x3)./h3; 
w4 = (P-x4)./h4;

xmid = (x1+x2+x3+x4)/4;
W = zeros(4,15);
W(1,:) = [double(subs(w1, P, x1 )),double(subs(w1, P, x2 )),double(subs(w1, P, x3 )),double(subs(w1, P, x4 )),double(subs(w1, P, xmid ))];
W(2,:) = [double(subs(w2, P, x1 )),double(subs(w2, P, x2 )),double(subs(w2, P, x3 )),double(subs(w2, P, x4 )),double(subs(w2, P, xmid ))];
W(3,:) = [double(subs(w3, P, x1 )),double(subs(w3, P, x2 )),double(subs(w3, P, x3 )),double(subs(w3, P, x4 )),double(subs(w3, P, xmid ))];
W(4,:) = [double(subs(w4, P, x1 )),double(subs(w4, P, x2 )),double(subs(w4, P, x3 )),double(subs(w4, P, x4 )),double(subs(w4, P, xmid ))];

% step 3: generating local MHFE matrices.

M = zeros(4,4);
for i = 1:4
    for j = 1:4
        M(i,j) = tdQuad(W(i,:), W(j,:), vol);
    end
end
condRT0 = cond(M);