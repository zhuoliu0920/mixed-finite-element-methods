function v = tdQuad(W1, W2, vol)
V = zeros(1,5);
for i=1:5
    V(i) = dot(W1((3*i-2):3*i), W2((3*i-2):3*i));
end
v = vol*(V(1)+V(2)+V(3)+V(4)+16*V(5))/20;
end

