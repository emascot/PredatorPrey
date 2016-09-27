data = importdata('data.txt');

x = 30/400:30/400:30;
y = 10/400:10/400:10;
[xx, yy] = meshgrid(x,y);
zz = griddata(data(:,1), data(:,2), data(:,3), xx, yy);

ax = gca;
surf(ax,xx,yy,zz,'EdgeColor','none');
ax.View = [0 90];
ax.XLabel.String = 'x (v_{prey}s)';
ax.YLabel.String = 'y (v_{prey}s)';
ax.ZLabel.String = 't (s)';
alpha(0.4);