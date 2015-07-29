clear;
n = 2;

if n > 1
seq(1) = 2;
row(1) = 1;
pos(1) = 2*n+2;
for i = 1:(2*n-1)
    for j = 2:(2*n+i)
        seq(2*n*(i-1)+i*(i-1)/2+j) = seq(2*n*(i-1)+i*(i-1)/2+j-1)+2;
        row(2*n*(i-1)+i*(i-1)/2+j) = row(2*n*(i-1)+i*(i-1)/2+j-1);
        if i == 1
            pos(2*n*(i-1)+i*(i-1)/2+j) = pos(2*n*(i-1)+i*(i-1)/2+j-1)+1;
        else
            pos(2*n*(i-1)+i*(i-1)/2+j) = pos(2*n*(i-1)+i*(i-1)/2+j-1);
        end
    end
    seq(2*n*i+i*(i+1)/2+1) = seq(2*n*i+i*(i+1)/2)+3;
    row(2*n*i+i*(i+1)/2+1) = row(2*n*i+i*(i+1)/2)+1;
    pos(2*n*i+i*(i+1)/2+1) = pos(2*n*i+i*(i+1)/2)+2;
end
end
