clear;
clc;

hold on
x=[0,1];
y=[0,1];
plot(x,y,'k');

x=[1,1];
y=[0,1];
plot(x,y,'k');

X=[0,0.5,1];
Y=[0,-0.5,0];

x=X(1:2);
y=Y(1:2);
plot(x,y,'k');

x=X(2:3);
y=Y(2:3);
plot(x,y,'k');

x=X(1:2:3);
y=Y(1:2:3);
plot(x,y,'k');

x=[X(2),(X(1)+X(3))/2];
y=[Y(2),(Y(1)+Y(3))/2];
plot(x,y,'k--');

hold off
axis equal
set(gca,'xlim',[min(X)-0.5 max(X)+0.5]);
set(gca,'ylim',[-1 1.5]);
set(gca,'Box','off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);