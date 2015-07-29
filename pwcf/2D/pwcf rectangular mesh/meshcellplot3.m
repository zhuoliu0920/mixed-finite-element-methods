clear;
clc;

x=[0,1/2,3/2];
y=[0,0,0];
plot(x,y,'k');
hold on

x=[0,1/2];
y=[1/2,1/2];
plot(x,y,'k');

x=[0,1/2,3/2];
y=[1,1,1];
plot(x,y,'k');

x=[0,0,0];
y=[0,1/2,1];
plot(x,y,'k');

x=[1/2,1/2,1/2];
y=[0,1/2,1];
plot(x,y,'k');

x=[3/2,3/2];
y=[0,1];
plot(x,y,'k');

x=[1/2,3/2];
y=[1,0];
plot(x,y,'k');

x=[0,1/2];
y=[1/2,0];
plot(x,y,'k');

x=[0,1/2];
y=[1/2,1];
plot(x,y,'k');

axis equal;
set(gca,'xlim',[0 3/2]);
set(gca,'ylim',[0 1]);
set(gca,'Box','off');
set(gca,'xtick',[])
set(gca,'ytick',[])