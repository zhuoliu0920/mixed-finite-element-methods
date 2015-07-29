clear;

hold on

X=[0,-1.25,1,0];
Y=[2,0.75,0.25,0];

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'--');

x=X(1:2);
y=Y(1:2);
plot(x,y,'--');

x=X(1:3:4);
y=Y(1:3:4);
plot(x,y,'k');

x=X(2:2:4);
y=Y(2:2:4);
plot(x,y,'k');

x=X(3:4);
y=Y(3:4);
plot(x,y,'k');

hold off
axis equal
axis off