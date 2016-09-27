% Import data
t = importdata('t.txt');
x = importdata('x.txt');
y = importdata('y.txt');

% Plot
ax = gca;
plot3(ax,x,y,t);
ax.View = [0 90];
ax.Title.String  = 'Predator Position';
ax.XLabel.String = 'x (v_{prey}s)';
ax.YLabel.String = 'y (v_{prey}s)';
ax.ZLabel.String = 't (s)';