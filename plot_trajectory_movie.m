% Import data
filename = 'trajectories.mp4';
t = importdata('t.txt');
x = importdata('x.txt');
y = importdata('y.txt');
N = size(x,2);
RGB = hsv(N);

% Initialize figure
fig = figure;
ax = gca;
fig.Position = [0 0 1440 800];

% Plot
hold(ax, 'on')
for i=1:N
    plot3(ax,x(:,i),y(:,i),t,'Color',RGB(i,:));
end
hold(ax, 'off')
ax.View = [0 90];
ax.Title.String  = 'Predator Position';
ax.XLabel.String = 'x (v_{prey}s)';
ax.YLabel.String = 'y (v_{prey}s)';
ax.ZLabel.String = 't (s)';

% Start recording
v = VideoWriter(filename,'MPEG-4');
v.FrameRate = 1000/30;
open(v);
for i=2:size(t);
    % Update
    ax.ZLim = [0 t(i)];
    frame = getframe(fig);
    writeVideo(v,frame);
end
close(v);