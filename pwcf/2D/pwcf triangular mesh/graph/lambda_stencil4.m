clear;
clc;

hold on

X=[-1,0,1];
Y=[0,sqrt(3),0];

x=X(1:2);
y=Y(1:2);
plot(x,y,'k');

x=X(2:3);
y=Y(2:3);
plot(x,y,'k');

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'k');

X=[-1,-0.5,0];
Y=[0,-sqrt(3)/2,0];

x=X(1:2);
y=Y(1:2);
plot(x,y,'k');

x=X(2:3);
y=Y(2:3);
plot(x,y,'k');

X=[1,0.5,0];
Y=[0,-sqrt(3)/2,0];

x=X(1:2);
y=Y(1:2);
plot(x,y,'k');

x=X(2:3);
y=Y(2:3);
plot(x,y,'k');

X=[-1.5,-1,-0.5];
Y=[-sqrt(3)/2,0,-sqrt(3)/2];

x=X(1:2);
y=Y(1:2);
plot(x,y,'k');

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'k');

hold off
axis equal
set(gca,'xlim',[-1.7 1.2]);
set(gca,'ylim',[-1.1 2]);
set(gca,'Box','off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);