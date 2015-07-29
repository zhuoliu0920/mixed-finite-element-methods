clear;
clc;

hold on

X=[0,0.3,0.75];
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

X=[0.75,1.05,1.5];
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
set(gca,'xlim',[-0.2 1.7]);
set(gca,'ylim',[-0.7 0.2]);
set(gca,'Box','off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);