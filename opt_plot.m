clc
clear

pos = [-1 -1 2 2];
rectangle('Position', pos, 'Curvature',[1 1], 'LineWidth', 0.5)
axis equal
%grid on
hold on

%% Heat
% plot(-0.2043, 0.978906,'ro', 'MarkerSize', 10)
% plot(-0.2043, 0.978906,'r.', 'MarkerSize', 14)
% plot(0.2043, -0.978906,'ro', 'MarkerSize', 10)
% plot(0.2043, -0.978906,'r.', 'MarkerSize', 14)
% plot(0.978906, -0.2043, 'ro', 'MarkerSize', 10)
% plot(0.978906, -0.2043, 'r.', 'MarkerSize', 14)
% plot(-0.978906, 0.2043, 'ro', 'MarkerSize', 10)
% plot(-0.978906, 0.2043, 'r.', 'MarkerSize', 14)

%% Wave
% plot(-0.0, 1.0,'bo', 'MarkerSize', 10)
% plot(-0.0, 1.0,'b.', 'MarkerSize', 14)
% plot(1.0, 0.0,'bo', 'MarkerSize', 10)
% plot(1.0, 0.0,'b.', 'MarkerSize', 14)
% plot(-1.0, 0.0, 'bo', 'MarkerSize', 10)
% plot(-1.0, 0.0,'b.', 'MarkerSize', 14)
% plot(-0.0, -1.0, 'bo', 'MarkerSize', 10)
% plot(-0.0, -1.0, 'b.', 'MarkerSize', 14)

%% Convection +
plot(-0.952449, 0.304699, 'mo', 'MarkerSize', 10)
plot(-0.952449, 0.304699, 'm.', 'MarkerSize', 14)
plot(0.952449, -0.304699, 'mo', 'MarkerSize', 10)
plot(0.952449, -0.304699, 'm.', 'MarkerSize', 14)

%% Convection - 
% plot(0.304699, -0.952449, 'mo', 'MarkerSize', 10)
% plot(0.304699, -0.952449, 'm.', 'MarkerSize', 14)
% plot(-0.304699, 0.952449, 'mo', 'MarkerSize', 10)
% plot(-0.304699, 0.952449, 'm.', 'MarkerSize', 14)

%% Plots
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
% For schrodinger:
%zlim([0 0.35]);
exportgraphics(ax,'convect_+opt.pdf','ContentType','vector')