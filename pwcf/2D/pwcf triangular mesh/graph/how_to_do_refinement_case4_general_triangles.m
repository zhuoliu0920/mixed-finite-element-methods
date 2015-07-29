clear;
clc;

hold on

theta1 = 80*pi/180;
theta2 = 75*pi/180;
t1 = tan(theta1); t2 = tan(theta2);
s1 = sin(theta1); s2 = sin(theta2);
c1 = cos(theta1); c2 = cos(theta2);
beta1 = pi*theta1/(theta1+theta2)/2;
beta2 = pi*theta2/(theta1+theta2)/2;

%{
ang=0:0.01:2*pi;
r = (1/t1+1/t2)/2; xp = (1/t2-1/t1)/2;
x = xp + r*cos(ang);
y = r*sin(ang);
plot(x,y,'r');
%}

%--------------------------------------------------------------------------

X=[-1/t1-1/t2,1/t2-1/t1,1/t1+1/t2];
Y=[-1,1,-1];

x=X(1:2);
y=Y(1:2);
plot(x,y,'k');

x=X(2:3);
y=Y(2:3);
plot(x,y,'k');

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'k');

%--------------------------------------------------------------------------

X=[-1/t1,0,1/t2];
Y=[0,-1,0];

x=X(1:2);
y=Y(1:2);
plot(x,y,'k');

x=X(2:3);
y=Y(2:3);
plot(x,y,'k');

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'k');

%--------------------------------------------------------------------------

a = sin(beta2)*(1/t1+1/t2);
b = sin(beta1)*(1/t1+1/t2);
ox = cos(beta1)*sin(beta2)*(1/t1+1/t2)-1/t1;
oy = -sin(beta1)*sin(beta2)*(1/t1+1/t2);

X=[-1/t1,ox,1/t2];
Y=[0,oy,0];

x=X(1:2);
y=Y(1:2);
plot(x,y,'--');

x=X(2:3);
y=Y(2:3);
plot(x,y,'--');

%--------------------------------------------------------------------------

X=[ox-a*sin(theta1-beta1)*sin(theta1),ox,ox+b*sin(theta2-beta2)*sin(theta2)];
Y=[oy-a*sin(theta1-beta1)*cos(theta1),oy,oy-b*sin(theta2-beta2)*cos(theta2)];

x=X(1:2);
y=Y(1:2);
plot(x,y,'--');

x=X(2:3);
y=Y(2:3);
plot(x,y,'--');

%--------------------------------------------------------------------------

x=[0,ox];
y=[-1,oy];
plot(x,y,'--');

%--------------------------------------------------------------------------

X=[ox-a*sin(theta1-beta1)*sin(theta1),-(1/t1+1/t2)/2,ox+b*sin(theta2-beta2)*sin(theta2)-(1/t1+1/t2)];
Y=[oy-a*sin(theta1-beta1)*cos(theta1),-1,oy-b*sin(theta2-beta2)*cos(theta2)];

x=X(1:2);
y=Y(1:2);
plot(x,y,'--');

x=X(2:3);
y=Y(2:3);
plot(x,y,'--');

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'--');

%--------------------------------------------------------------------------

X=[ox-a*sin(theta1-beta1)*sin(theta1)+(1/t1+1/t2),(1/t1+1/t2)/2,ox+b*sin(theta2-beta2)*sin(theta2)];
Y=[oy-a*sin(theta1-beta1)*cos(theta1),-1,oy-b*sin(theta2-beta2)*cos(theta2)];

x=X(1:2);
y=Y(1:2);
plot(x,y,'--');

x=X(2:3);
y=Y(2:3);
plot(x,y,'--');

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'--');

%--------------------------------------------------------------------------

hold off
axis equal
set(gca,'xlim',[-(1/t1+1/t2)-0.2 (1/t1+1/t2)+0.2]);
set(gca,'ylim',[-1.2 1.2]);
set(gca,'Box','off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);

% check if there is any abtuse triangle

v1 = [X(1)-X(2),Y(1)-Y(2)];
v2 = [X(3)-X(1),Y(3)-Y(1)];
v3 = [X(3)-X(2),Y(3)-Y(2)];
v4 = [1/t1,-1];
v5 = [-1/t2,-1];
angle1 = acos(-dot(v1,v2)/(norm(v1)*norm(v2)))*180/pi;
angle2 = acos(dot(v1,v3)/(norm(v1)*norm(v3)))*180/pi;
angle3 = acos(dot(v3,v2)/(norm(v3)*norm(v2)))*180/pi;
angle4 = acos(-dot(v4,v2)/(norm(v4)*norm(v2)))*180/pi;
angle5 = acos(dot(v5,v2)/(norm(v5)*norm(v2)))*180/pi;