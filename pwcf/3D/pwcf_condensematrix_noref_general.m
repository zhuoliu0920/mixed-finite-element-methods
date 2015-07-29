clear;
% input coordinates of four vertexes.

h = 1;

%x1 = [0,0,1]; x4 = [0,0,0]; x2 = [0,1,0]; x3 = [1,0,0];

%x1 = h.*[sqrt(3)/3,0,2*sqrt(6)/3]; x2 = h.*[0,-1,0]; x3 = h.*[0,1,0]; x4 = h.*[sqrt(3),0,0];

%x1 = [1,1,1]; x4 = [0,0,0]; x2 = [1,-1,0]; x3 = [2,2,0]; 

x1 = [0,0,100]; x2 = [-1,0,0]; x3 = [1,0,0]; x4 = [0,1,0]; 

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
B = -[s1,s2,s3,s4];
C = diag(-B);
invM = inv(M);
b = B*inv(M)*B';
A = inv(M)*B';
D = (1/b).*(A*A');
E = inv(M)-D;
S = C*E*C';

h11 = B*inv(M)*B'; h12 = B*inv(M)*C'; h21 = C*inv(M)*B'; h22 = C*inv(M)*C';
S2 = [h11 h12; h21 h22];

S3 = h21*h12/h11;

alpha = dot(n2, cross(n1,n4))/norm(cross(n1,n4));
beta = dot(n3, cross(n1,n4))/norm(cross(n1,n4));
j11 = -4*alpha*beta*s2*s3/vol;
j12 = 1/vol.*[0, 2*alpha*beta*s2*s3, 2*alpha*beta*s2*s3, 0];
j21 = j12';
j22 = S; 
j22(2,2) = j22(2,2)+alpha*alpha*s2*s2/vol;
j22(2,3) = j22(2,3)-alpha*beta*s2*s3/vol;
j22(3,2) = j22(2,3);
j22(3,3) = j22(3,3)+beta*beta*s3*s3/vol;
S4 = [j11 j12; j21 j22];









%matrix2latex(M, '3Mp.tex');
%matrix2latex(B, '3Bp.tex');
%matrix2latex(C, '3Cp.tex');
%matrix2latex(invM, '3invMp.tex');
%matrix2latex(S, '3Sp.tex');
%matrix2latex(S2, '3S2p.tex');
%matrix2latex(S3, '3S3p.tex');





%m = cross(n1,n4);

%S1 = zeros(4,4);
%S1(1,1) = dot(n1,n1)*s1*s1; S1(1,2) = dot(n1,n2)*s1*s2; S1(1,3) = dot(n1,n3)*s1*s3; S1(1,4) = dot(n1,n4)*s1*s4;
%S1(2,2) = dot(n2,n2)*s2*s2; S1(2,3) = dot(n2,n3)*s2*s3; S1(2,4) = dot(n2,n4)*s2*s4;
%S1(3,3) = dot(n3,n3)*s3*s3; S1(3,4) = dot(n3,n4)*s3*s4;
%S1(4,4) = dot(n4,n4)*s4*s4;
%S1(2,1) = S1(1,2); S1(3,1) = S1(1,3); S1(3,2) = S1(2,3); S1(4,1) = S1(1,4); S1(4,2) = S1(2,4); S1(4,3) = S1(3,4);
%S1 = 1/vol.*S1;