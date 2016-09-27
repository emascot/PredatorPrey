t = importdata('t.txt');
x = importdata('x.txt');
y = importdata('y.txt');

ax = gca;
plot3(ax,x,y,t);
ax.View = [0 90];
ax.XLabel.String = 'x (v_{prey}s)';
ax.YLabel.String = 'y (v_{prey}s)';
ax.ZLabel.String = 't (s)';