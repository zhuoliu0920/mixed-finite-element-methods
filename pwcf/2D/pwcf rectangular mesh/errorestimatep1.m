function [e,preal] = errorestimatep1(p,n,stage)
h = 1/n;
preal = zeros(n*n+3*stage,1);
for i = 1:(n-1)
    for j = 1:n
        preal(j+(i-1)*n) = 2-(i+j-1)/n;
    end
end
if stage == 0
    i = n;
    for j = 1:n
        preal(j+(i-1)*n) = 2-(i+j-1)/n;
    end
else
    i = n;
    for j = 1:(n-1)
        preal(j+(i-1)*n) = 2-(i+j-1)/n;
    end
    for s = 1:stage
        preal(n*n+3*(s-1)) = 3/(n*2^s);
        preal(n*n+3*(s-1)+1) = 2/(n*2^s);
        preal(n*n+3*(s-1)+2) = 2/(n*2^s);
    end
    preal(n*n+3*stage) = 1/(n*2^stage);
end
e = max(abs(p-preal));