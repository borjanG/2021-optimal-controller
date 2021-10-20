clc
clear

pos = [-1 -1 2 2];
rectangle('Position', pos, 'Curvature',[1 1], 'LineWidth', 0.5)
axis equal
%grid on
hold on

%% Heat
% plot(-0.25798309, 0.96614944,'ro', 'MarkerSize', 18)
% plot(-0.25798309, 0.96614944,'r.', 'MarkerSize', 24)
% plot(0.25798309, -0.96614944,'ro', 'MarkerSize', 18)
% plot(0.25798309, -0.96614944,'r.', 'MarkerSize', 24)
% plot(0.96614944, -0.25798309, 'ro', 'MarkerSize', 18)
% plot(0.96614944, -0.25798309, 'r.', 'MarkerSize', 24)
% plot(-0.96614944, 0.25798309, 'ro', 'MarkerSize', 18)
% plot(-0.96614944, 0.25798309, 'r.', 'MarkerSize', 24)

%% Wave
% plot(-0.25798309, 0.96614944,'bo', 'MarkerSize', 18)
% plot(-0.25798309, 0.96614944,'b.', 'MarkerSize', 24)
% plot(0.25798309, -0.96614944,'bo', 'MarkerSize', 18)
% plot(0.25798309, -0.96614944,'b.', 'MarkerSize', 24)
% plot(0.96614944, -0.25798309, 'bo', 'MarkerSize', 18)
% plot(0.96614944, -0.25798309, 'b.', 'MarkerSize', 24)
% plot(-0.96614944, 0.25798309, 'bo', 'MarkerSize', 18)
% plot(-0.96614944, 0.25798309, 'b.', 'MarkerSize', 24)

%% Convection +-
% plot(-0.9549099, 0.29689575, 'mo', 'MarkerSize', 18)
% plot(-0.9549099, 0.29689575, 'm.', 'MarkerSize', 24)
% plot(0.9549099, -0.29689575, 'mo', 'MarkerSize', 18)
% plot(0.9549099, -0.29689575, 'm.', 'MarkerSize', 24)

%% Convection -+ 
plot(0.29689575, -0.9549099, 'mo', 'MarkerSize', 18)
plot(0.29689575, -0.9549099, 'm.', 'MarkerSize', 24)
plot(-0.29689575, 0.9549099, 'mo', 'MarkerSize', 18)
plot(-0.29689575, 0.9549099, 'm.', 'MarkerSize', 24)

%% Plots
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
set(gca,'XMinorTick','on','YMinorTick','on')
grid minor
% For schrodinger:
%zlim([0 0.35]);
exportgraphics(ax,'heat_opt.pdf','ContentType','vector')