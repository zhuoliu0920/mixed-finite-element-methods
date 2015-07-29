clear;
n = 16;
[u,p,lambda, M, B, C, H,Global] = pwcf_tri_noref(n);
preal = pwcf_real_sol_noref(n);
error = 0;
for i = 1:n*n
    error = error + sqrt(3)/(4*n*n)*(preal(i)-p(i))^(2);
end
error = sqrt(error);