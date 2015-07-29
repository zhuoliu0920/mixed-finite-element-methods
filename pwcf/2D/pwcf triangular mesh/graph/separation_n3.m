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

for i = 1:3
    x(1) = -i/6;
    x(2) = i/6;
    y(1) = sqrt(3)-i*sqrt(3)/6;
    y(2) = y(1);
    plot(x,y,'k');

    x(1) = -i/6;
    x(2) = 0.5-i/3;
    y(1) = sqrt(3)-i*sqrt(3)/6;
    y(2) = sqrt(3)/2;
    plot(x,y,'k');
    
    x(1) = i/6;
    x(2) = -0.5+i/3;
    y(1) = sqrt(3)-i*sqrt(3)/6;
    y(2) = sqrt(3)/2;
    plot(x,y,'k');
end

for i = 0:6
    x(1) = -0.5-i/12;
    x(2) = 0.5+i/12;
    y(1) = sqrt(3)/2-i*sqrt(3)/12;
    y(2) = y(1);
    plot(x,y,'k');
    
    x(1) = -0.5-i/12;
    x(2) = -i/6;
    y(1) = sqrt(3)/2-i*sqrt(3)/12;
    y(2) = 0;
    plot(x,y,'k');
    
    x(1) = 0.5+i/12;
    x(2) = i/6;
    y(1) = sqrt(3)/2-i*sqrt(3)/12;
    y(2) = 0;
    plot(x,y,'k');
end

for i = 1:5
    x(1) = -0.5+i/6;
    x(2) = i/6;
    y(1) = sqrt(3)/2;
    y(2) = 0;
    plot(x,y,'k');
    
    x(1) = -0.5+i/6;
    x(2) = -1+i/6;
    y(1) = sqrt(3)/2;
    y(2) = 0;
    plot(x,y,'k');
end

hold off
axis equal
set(gca,'xlim',[-1.2 1.2]);
set(gca,'ylim',[-2 2]);
set(gca,'Box','off');
set(gca,'xtick',[]);
set(gca,'ytick',[]);