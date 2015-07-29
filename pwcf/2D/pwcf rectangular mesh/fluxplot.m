clear all
[X,Y] = meshgrid(1/120:1/60:1);
Z = sin(pi.*X).*sin(pi.*Y);
[DX,DY] = gradient(Z);
quiver(X,Y,DX,DY);
grid on;
set(gca,'xtick',0:1/30:1);
set(gca,'ytick',0:1/30:1);
set(gca,'xlim',[0 1]);
set(gca,'ylim',[0 1]);