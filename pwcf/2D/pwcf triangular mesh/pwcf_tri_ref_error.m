clear;
n = 8;
[u,p,lambda, M, B, C, H,Global] = pwcf_tri_ref(n);
preal = pwcf_real_sol(n);
error = 0;
for i = 1:n*n
    error = error + sqrt(3)/(4*n*n)*(preal(i)-p(i))^(2);
end
for i = (n*n+1):(13*n*n)
    error = error + sqrt(3)/(16*n*n)*(preal(i)-p(i))^(2);
end
error = sqrt(error);