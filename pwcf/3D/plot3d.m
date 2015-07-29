clear;
clc;

X=[0,2,0.5,1.5,1.25];
Y=[0,1,2,1,1.5];
Z=[0,0,0,1,0];

x=X(1:2);
y=Y(1:2);
z=Z(1:2);
plot3(x,y,z,'k');
hold on

x=X(1:2:3);
y=Y(1:2:3);
z=Z(1:2:3);
plot3(x,y,z,'k');

x=X(1:3:4);
y=Y(1:3:4);
z=Z(1:3:4);
plot3(x,y,z,'k');

x=X(2:3);
y=Y(2:3);
z=Z(2:3);
plot3(x,y,z,'k--');

x=X(2:2:4);
y=Y(2:2:4);
z=Z(2:2:4);
plot3(x,y,z,'k');

x=X(3:4);
y=Y(3:4);
z=Z(3:4);
plot3(x,y,z,'k');

x=X(1:4:5);
y=Y(1:4:5);
z=Z(1:4:5);
plot3(x,y,z,'b:','LineWidth',2);

x=X(4:5);
y=Y(4:5);
z=Z(4:5);
plot3(x,y,z,'b:','LineWidth',2);

X=[X;Y;Z];

text(X(1,1),X(2,1)-0.05,X(3,1),'v1','FontSize',12)
text(X(1,2)+0.02,X(2,2),X(3,2),'v2','FontSize',12) 
text(X(1,3)-0.1,X(2,3),X(3,3),'v3','FontSize',12) 
text(X(1,4)+0.02,X(2,4),X(3,4),'v4','FontSize',12)

text((X(1,1)+X(1,3)+X(1,4))/3-0.2,(X(2,1)+X(2,3)+X(2,4))/3,(X(3,1)+X(3,3)+X(3,4))/3,'e_{1}','FontSize',15,'FontWeight','bold') 
text((X(1,1)+X(1,2)+X(1,4))/3,(X(2,1)+X(2,2)+X(2,4))/3,(X(3,1)+X(3,2)+X(3,4))/3,'e_{2}','FontSize',15,'FontWeight','bold') 

hold off
axis equal
axis off

print -depsc tetrahedron.eps