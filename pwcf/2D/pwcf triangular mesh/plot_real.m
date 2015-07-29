%X = linspace(-1,1,50);
%Y = linspace(0,sqrt(3),50);
%[x y] = meshgrid(X,Y);
x = 0; y = sqrt(3)/2;
z = -y.*(sqrt(3).*x+y-sqrt(3)).*(sqrt(3).*x-y+sqrt(3));
%surf(x,y,z)