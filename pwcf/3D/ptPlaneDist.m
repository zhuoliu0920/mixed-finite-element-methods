function h = ptPlaneDist(P1,P2,P3,P4 )
normal = cross(P2-P3, P2-P4);
syms x y z;
P = [x,y,z];
planefunction = dot(normal, P-P2);
value = double( subs(planefunction, P, P1 ) );
h = abs(value/norm(normal));
end

