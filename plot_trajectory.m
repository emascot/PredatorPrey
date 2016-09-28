% Import data
t = importdata('t.txt');
x = importdata('x.txt');
y = importdata('y.txt');
N = size(x,2);
RGB = hsv(N);

% Plot
ax = gca;
hold(ax, 'on')
for i=1:1000
    plot3(ax,x(:,i),y(:,i),t,'Color',RGB(i,:));
end
hold(ax, 'off')

ax.View = [0 90];
ax.Title.String  = 'Predator Position';
ax.XLabel.String = 'x (v_{prey}s)';
ax.YLabel.String = 'y (v_{prey}s)';
ax.ZLabel.String = 't (s)';
ax.YLim = [0.01 10];