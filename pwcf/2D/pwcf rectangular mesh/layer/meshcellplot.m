clear;
clc;
hold on
x=[0:0.5:10];
y=[0:0.5:10];
M=meshgrid(x,y);
plot(x,M,'k');

for i = -12:12
x=[0+i,10+i];
y=[0,10];
plot(1-x,y,'k');

x=[0+i,10+i];
y=[0,10];
plot(x,y,'k');

x=[0.5+i,2+i];
y=[0,1.5];
plot(x,y,'k');

x=[0.5+i,2+i];
y=[0,1.5];
plot(1-x,y,'k');

x=[0,10];
y=[0.25,0.25];
plot(x,y,'k');

x=[0,10];
y=[0.75,0.75];
plot(x,y,'k');

x=[0,10];
y=[1.25,1.25];
plot(x,y,'k');
end
hold off

axis equal;
set(gca,'xlim',[0 10]);
set(gca,'ylim',[0 3]);
set(gca,'Box','off');
set(gca,'xtick',[])
set(gca,'ytick',[])