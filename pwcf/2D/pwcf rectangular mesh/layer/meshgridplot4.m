clear;
clc;

hold on
x=[0,1];
y=[0,1];
plot(x,y,'k');

x=[1.5,1];
y=[0,1];
plot(x,y,'k');

X=[0,0.375,1.5];
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

hold off
axis equal
set(gca,'xlim',[-0.2 1.7]);
set(gca,'ylim',[-1.2 1.2]);
set(gca,'Box','off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);