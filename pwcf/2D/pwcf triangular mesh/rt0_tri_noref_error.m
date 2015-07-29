n = 16;
[u,p,lambda, M, B, C, H,Global] = rt0_tri_noref(n);
preal = rt0_real_sol_noref(n);
error = 0;
for i = 1:n*n
    error = error + sqrt(3)/(n*n)*(preal(i)-p(i))^(2);
end
error = sqrt(error);