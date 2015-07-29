function p = pwcf_real_sol(n)

h = sqrt(3)/2/n;
b = 1/n;
row_ele(1) = 1;
total_ele(1) = 1;
for i = 2:n
    row_ele(i) = 2*i-1;
    total_ele(i) = total_ele(i-1)+row_ele(i);
end
for i = (n+1):3*n
    row_ele(i) = 4*n-1+2*(i-n);
    total_ele(i) = total_ele(i-1)+row_ele(i);
end

p = zeros(13*n*n,1);

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

i = n+1;
for j = 1:2*n
    x = -0.5+b/4+(j-1)*(b/2);
    fun = @(y) (y.*y.*y-2*sqrt(3).*y.*y+3.*y-3*x*x.*y);
    p(j+total_ele(i-1)) = quad(fun,sqrt(3)-n*h-(i-n)*h/2,sqrt(3)-n*h-(i-n-1)*h/2)/(h/2);
end
for j = (2*n+1):(4*n+1)
    x = -0.5+(j-2*n-1)*(b/2);
    fun = @(y) (y.*y.*y-2*sqrt(3).*y.*y+3.*y-3*x*x.*y);
    p(j+total_ele(i-1)) = quad(fun,sqrt(3)-n*h-(i-n)*h/2,sqrt(3)-n*h-(i-n-1)*h/2)/(h/2);
end
    
for i = (n+2):3*n
    for j =1:row_ele(i)
        x = (j-(row_ele(i)+1)/2)*(b/4);
        fun = @(y) (y.*y.*y-2*sqrt(3).*y.*y+3.*y-3*x*x.*y);
        p(j+total_ele(i-1)) = quad(fun,sqrt(3)-n*h-(i-n)*h/2,sqrt(3)-n*h-(i-n-1)*h/2)/(h/2);
    end
end
    


end

