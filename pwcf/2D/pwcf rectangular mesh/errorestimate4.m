function [e,uh] = errorestimate4(u,n,stage)
% g = x1^2 + x2^2, f = 0;
uh = zeros(4,1);
e = 0;
for j = 1:(n-1)
    for i = 1:n
        uh(1) = -u(2+4*(i+n*(j-1)-1));
        uh(2) = -u(1+4*(i+n*(j-1)-1));
        uh(3) = -u(4+4*(i+n*(j-1)-1));
        uh(4) = -u(3+4*(i+n*(j-1)-1));
        
        fun = @(x1,x2) (2*x1+uh(1)).^2+(2*x2+uh(2)).^2;
        x1min = 1-j/n; x1max = 1-(j-1)/n; x2min = 1-i/n; x2max = @(x1) -x1+2-(i+j-1)/n;
        e = e + quad2d(fun,x1min,x1max,x2min,x2max);
        fun = @(x1,x2) (2*x1+uh(3)).^2+(2*x2+uh(4)).^2;
        x2min = @(x1) -x1+2-(i+j-1)/n; x2max = 1-(i-1)/n;
        e = e + quad2d(fun,x1min,x1max,x2min,x2max);
    end
end
j = n;
if stage == 0
   for i = 1:n
       uh(1) = -u(2+4*(i+n*(j-1)-1));
       uh(2) = -u(1+4*(i+n*(j-1)-1));
       uh(3) = -u(4+4*(i+n*(j-1)-1));
       uh(4) = -u(3+4*(i+n*(j-1)-1));
        
       fun = @(x1,x2) (2*x1+uh(1)).^2+(2*x2+uh(2)).^2;
       x1min = 1-j/n; x1max = 1-(j-1)/n; x2min = 1-i/n; x2max = @(x1) -x1+2-(i+j-1)/n;
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
       fun = @(x1,x2) (2*x1+uh(3)).^2+(2*x2+uh(4)).^2;
       x2min = @(x1) -x1+2-(i+j-1)/n; x2max = 1-(i-1)/n;
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
   end

else
   for i = 1:(n-1)
       uh(1) = -u(2+4*(i+n*(j-1)-1));
       uh(2) = -u(1+4*(i+n*(j-1)-1));
       uh(3) = -u(4+4*(i+n*(j-1)-1));
       uh(4) = -u(3+4*(i+n*(j-1)-1));
        
       fun = @(x1,x2) (2*x1+uh(1)).^2+(2*x2+uh(2)).^2;
       x1min = 1-j/n; x1max = 1-(j-1)/n; x2min = 1-i/n; x2max = @(x1) -x1+2-(i+j-1)/n;
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
       fun = @(x1,x2) (2*x1+uh(3)).^2+(2*x2+uh(4)).^2;
       x2min = @(x1) -x1+2-(i+j-1)/n; x2max = 1-(i-1)/n;
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
   end
   
   for s = 1:stage % local refinement
       
       i = 1;
       uh(1) = -u(2+4*(n*n+3*(s-1)+i-1));
       uh(2) = -u(1+4*(n*n+3*(s-1)+i-1));
       uh(3) = -u(4+4*(n*n+3*(s-1)+i-1));
       uh(4) = -u(3+4*(n*n+3*(s-1)+i-1));
       fun = @(x1,x2) (2*x1+uh(1)).^2+(2*x2+uh(4)).^2;
       x1min = 1/(n*2^s); x1max = 1/(n*2^(s-1)); x2min = @(x1) x1; x2max = 1/(n*2^(s-1));
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
       fun = @(x1,x2) (2*x1+uh(3)).^2+(2*x2+uh(2)).^2;
       x2min = 1/(n*2^s); x2max = @(x1) x1;
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
       
       i = 2;
       uh(1) = -u(2+4*(n*n+3*(s-1)+i-1));
       uh(2) = -u(1+4*(n*n+3*(s-1)+i-1));
       uh(3) = -u(4+4*(n*n+3*(s-1)+i-1));
       uh(4) = -u(3+4*(n*n+3*(s-1)+i-1));
       fun = @(x1,x2) (2*x1+uh(1)).^2+(2*x2+uh(2)).^2;
       x1min = 1/(n*2^s); x1max = 1/(n*2^(s-1)); x2min = (2-i)/(n*2^s); x2max = @(x1) -x1+x1max+x2min;
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
       fun = @(x1,x2) (2*x1+uh(3)).^2+(2*x2+uh(4)).^2;
       x2min = @(x1) -x1+x1max+x2min; x2max = 2/(i*(n*2^s));
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);      
       
       
       i = 3;
       uh(1) = -u(2+4*(n*n+3*(s-1)+i-1));
       uh(2) = -u(1+4*(n*n+3*(s-1)+i-1));
       uh(3) = -u(4+4*(n*n+3*(s-1)+i-1));
       uh(4) = -u(3+4*(n*n+3*(s-1)+i-1));
       
       fun = @(x1,x2) (2*x1+uh(1)).^2+(2*x2+uh(2)).^2;
       x1min = 0; x1max = 1/(n*2^s); x2min = 1/(n*2^s); x2max = @(x1) -x1+x1max+x2min;
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
       fun = @(x1,x2) (2*x1+uh(3)).^2+(2*x2+uh(4)).^2;
       x2min = @(x1) -x1+x1max+x2min; x2max = 1/(n*2^(s-1));
       e = e + quad2d(fun,x1min,x1max,x2min,x2max);
   end
   % the last mesh cell
   uh(1) = -u(2+4*(3*stage+n*n-1));
   uh(2) = -u(1+4*(3*stage+n*n-1));
   uh(3) = -u(4+4*(3*stage+n*n-1));
   uh(4) = -u(3+4*(3*stage+n*n-1));
        
   fun = @(x1,x2) (2*x1+uh(1)).^2+(2*x2+uh(4)).^2;
   x1min = 0; x1max = 1/(n*2^stage); x2min = @(x1) x1; x2max = 1/(n*2^stage);
   e = e + quad2d(fun,x1min,x1max,x2min,x2max);
   fun = @(x1,x2) (2*x1+uh(3)).^2+(2*x2+uh(2)).^2;
   x2min = 0; x2max = @(x1) x1;
   e = e + quad2d(fun,x1min,x1max,x2min,x2max);
end
e = sqrt(e);