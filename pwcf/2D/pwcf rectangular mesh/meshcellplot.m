clear;
clc;

x=[0,1/2,1];
y=[0,1/2,1];
M=meshgrid(x,y);
plot(x,M,'k');
hold on
plot(M,y,'k');

x=[0,1/4,1/2];
y=[0,1/4,1/2];
M=meshgrid(x,y);
plot(x,M,'k');
hold on
plot(M,y,'k');

x=[0,1/8,1/4];
y=[0,1/8,1/4];
M=meshgrid(x,y);
plot(x,M,'k');
hold on
plot(M,y,'k');

x=[0,1];
y=[0,1];
plot(1-x,y,'k');
plot(3/2-x,y,'k');
plot(1/2-x,y,'k');
plot(1/4-x,y,'k');

x=[0,1/2];
y=[0,1/2];
plot(x,y,'k');

axis equal;
set(gca,'xlim',[0 1]);
set(gca,'ylim',[0 1]);
set(gca,'Box','off');
set(gca,'xtick',[])
set(gca,'ytick',[])