function [e,preal] = errorestimatep(p,n,stage)
h = 1/n;
preal = zeros(n*n+3*stage,1);
for i = 1:(n-1)
    for j = 1:n
        preal(j+(i-1)*n) = (quad(@(x1) sqrt(2)*sin(pi.*x1).*sin(pi.*(2-(i+j-1)/n-x1)),1-i/n,1-(i-1)/n))/(sqrt(2)*h);
    end
end
if stage == 0
    i = n;
    for j = 1:n
        preal(j+(i-1)*n) = (quad(@(x1) sqrt(2)*sin(pi.*x1).*sin(pi.*(2-(i+j-1)/n-x1)),1-i/n,1-(i-1)/n))/(sqrt(2)*h);
    end
else
    i = n;
    for j = 1:(n-1)
        preal(j+(i-1)*n) = (quad(@(x1) sqrt(2)*sin(pi.*x1).*sin(pi.*(2-(i+j-1)/n-x1)),1-i/n,1-(i-1)/n))/(sqrt(2)*h);
    end
    for s = 1:stage
        preal(n*n+3*(s-1)) = (quad(@(x1) sqrt(2)*sin(pi.*x1).*sin(pi.*(3/(n.*2.^s)-x1)),1/(n*2^s),1/(n*2^(s-1))))/(sqrt(2)*h/2^s);
        preal(n*n+3*(s-1)+1) = (quad(@(x1) sqrt(2)*sin(pi.*x1).*sin(pi.*(2/(n.*2.^s)-x1)),1/(n*2^s),1/(n*2^(s-1))))/(sqrt(2)*h/2^s);
        preal(n*n+3*(s-1)+2) = (quad(@(x1) sqrt(2)*sin(pi.*x1).*sin(pi.*(2/(n.*2.^s)-x1)),0,1/(n*2^s)))/(sqrt(2)*h/2^s);
    end
    preal(n*n+3*stage) = (quad(@(x1) sqrt(2)*sin(pi.*x1).*sin(pi.*(1/(n.*2.^stage)-x1)),0,1/(n*2^stage)))/(sqrt(2)*h/2^stage);
end
e = max(abs(p-preal));