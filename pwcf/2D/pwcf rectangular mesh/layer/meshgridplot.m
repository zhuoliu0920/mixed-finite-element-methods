clear;
clc;

hold on

X=[0,0.5,1.5];
Y=[0,1,0];

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
set(gca,'xlim',[min(X)-0.3 max(X)+0.1]);
set(gca,'ylim',[min(Y)-0.3 max(Y)+0.1]);
set(gca,'Box','off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);