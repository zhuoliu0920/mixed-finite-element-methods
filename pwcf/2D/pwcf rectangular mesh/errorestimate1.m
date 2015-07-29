function [e,uh] = errorestimate1(u,n,stage)
% g = x1 + x2, f = 0;
h = 1/n;
uh = zeros(4,1);
e = 0;
for j = 1:(n-1)
    for i = 1:n
        uh(1) = -u(2+4*(i+n*(j-1)-1));
        uh(2) = -u(1+4*(i+n*(j-1)-1));
        uh(3) = -u(4+4*(i+n*(j-1)-1));
        uh(4) = -u(3+4*(i+n*(j-1)-1));
        e = e + h*h/2*((1+uh(1)).^2+(1+uh(2)).^2+(1+uh(3)).^2+(1+uh(4)).^2);
    end
end
j = n;
if stage == 0
   for i = 1:n
       uh(1) = -u(2+4*(i+n*(j-1)-1));
       uh(2) = -u(1+4*(i+n*(j-1)-1));
       uh(3) = -u(4+4*(i+n*(j-1)-1));
       uh(4) = -u(3+4*(i+n*(j-1)-1));
        
       e = e + h*h/2*((1+uh(1)).^2+(1+uh(2)).^2+(1+uh(3)).^2+(1+uh(4)).^2);
   end

else
   for i = 1:(n-1)
       uh(1) = -u(2+4*(i+n*(j-1)-1));
       uh(2) = -u(1+4*(i+n*(j-1)-1));
       uh(3) = -u(4+4*(i+n*(j-1)-1));
       uh(4) = -u(3+4*(i+n*(j-1)-1));
        
       e = e + h*h/2*((1+uh(1)).^2+(1+uh(2)).^2+(1+uh(3)).^2+(1+uh(4)).^2);
   end  
   for s = 1:stage % local refinement
       for i = 1:3
       uh(1) = -u(2+4*(n*n+3*(s-1)+i-1));
       uh(2) = -u(1+4*(n*n+3*(s-1)+i-1));
       uh(3) = -u(4+4*(n*n+3*(s-1)+i-1));
       uh(4) = -u(3+4*(n*n+3*(s-1)+i-1));
       
       e = e + h*h/(2*4^s)*((1+uh(1)).^2+(1+uh(2)).^2+(1+uh(3)).^2+(1+uh(4)).^2);
       end
   end
   % the last mesh cell
   uh(1) = -u(2+4*(3*stage+n*n-1));
   uh(2) = -u(1+4*(3*stage+n*n-1));
   uh(3) = -u(4+4*(3*stage+n*n-1));
   uh(4) = -u(3+4*(3*stage+n*n-1));
        
   e = e + h*h/(2*4^stage)*((1+uh(1)).^2+(1+uh(2)).^2+(1+uh(3)).^2+(1+uh(4)).^2);
end
e = sqrt(e);