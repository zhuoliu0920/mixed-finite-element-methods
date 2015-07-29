clear;
clc;
hold on
x=[0:0.5:2];
y=[0:0.5:2];
M=meshgrid(x,y);
plot(x,M,'k');

for i = -12:12

x=[0+i/2,10+i/2];
y=[0,10];
plot(x,y,'k');

x=[i/2,i/2];
y=[0,10];
plot(x,y,'k');


end
hold off

axis equal;
set(gca,'xlim',[0 2]);
set(gca,'ylim',[0 2]);
set(gca,'Box','off');
set(gca,'xtick',[])
set(gca,'ytick',[])