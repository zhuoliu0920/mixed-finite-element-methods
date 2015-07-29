function p = pwcf_real_sol_noref(n)

h = sqrt(3)/n;
b = 2/n;
row_ele(1) = 1;
total_ele(1) = 1;
for i = 2:n
    row_ele(i) = 2*i-1;
    total_ele(i) = total_ele(i-1)+row_ele(i);
end

p = zeros(n*n,1);

x = 0;
fun = @(y) (y.*y.*y-2*sqrt(3).*y.*y+3.*y-3*x*x.*y);
p(1) = quad(fun,sqrt(3)-h,sqrt(3))/h;

for i = 2:n
    for j =1:row_ele(i)
        x = (j-(row_ele(i)+1)/2)*(b/2);
        fun = @(y) (y.*y.*y-2*sqrt(3).*y.*y+3.*y-3*x*x.*y);
        p(j+total_ele(i-1)) = quad(fun,sqrt(3)-i*h,sqrt(3)-(i-1)*h)/h;
    end
end

end

