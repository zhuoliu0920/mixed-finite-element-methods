clear;
clc;
x=[0,1];
y=[0,1];
M=meshgrid(x,y);
plot(x,M,'k');
hold on
plot(M,y,'k');
plot(1-x,y,'k')
X = [0;1/2;1/2;1];
Y = [1/2;0;1;1/2];
DX = [-1/10;0;0;-1/10];
DY = [0;-1/10;-1/10;0];
quiver(X,Y,DX,DY,0.2);
axis square
set(gca,'xlim',[-0.3 1.1]);
set(gca,'ylim',[-0.3 1.1]);
set(gca,'Box','off')