o = pi/6;
a = sec(o); b = 0; c = 0;
em1 = [a*cot(o)/2 + b, c];
em2 = [b, a/2 + c];
em3 = [a*cot(o)/2 + b, a/2 + c];
n1 = [0;-1]; n2 = [-1;0]; n3=[sin(o);cos(o)];

result = [em1*n1;em2*n2;em3*n3];