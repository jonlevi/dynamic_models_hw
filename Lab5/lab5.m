%% Lab 5
% Part 1
% Pre-bifurcation
mu = 0.5;

dx_dt = @(x,y) mu-x.^2;
dy_dt = @(x,y) -y;
figure();
[x, y] = meshgrid(-2:0.02:2, -2:0.02:2);

dx = dx_dt(x, y);
dy = dy_dt(x, y);

streamslice(x, y, dx, dy);
xlabel('x'), ylabel('y')
title('\mu=0.5')
%axis('tight', 'equal')
hold on;

%% 
% Bifurcation
mu = 0;

dx_dt = @(x,y) mu-x.^2;
dy_dt = @(x,y) -y;
figure();
[x, y] = meshgrid(-2:0.02:2, -2:0.02:2);

dx = dx_dt(x, y);
dy = dy_dt(x, y);

streamslice(x, y, dx, dy);
xlabel('x'), ylabel('y')
%axis('tight', 'equal')
title('\mu=0')
hold on;
%% 
% Post- Bifurcation
mu = -0.5;

dx_dt = @(x,y) mu-x.^2;
dy_dt = @(x,y) -y;
figure();
[x, y] = meshgrid(-8:0.02:8, -8:0.02:8);

dx = dx_dt(x, y);
dy = dy_dt(x, y);

streamslice(x, y, dx, dy);
xlabel('x'), ylabel('y')
%axis('tight', 'equal')
title('\mu=-0.5')
hold on;
%% 
% Post- Bifurcation time series
mu = -0.5;

dx_dt = @(x,y) mu-x.^2;
dy_dt = @(x,y) -y;
figure();
options = odeset('RelTol',1e-8);
[t,ysolution] = ode45(@(t,y) [mu-(y(1).^2); -y(2)],[0 4.25],[10 -3],options);
[t2,ysolution2] = ode45(@(t,y) [mu-(y(1).^2); -y(2)],[0 10],[25 3],options);
subplot(1,2,1);
plot(t, ysolution(:,1), 'b');
ylim([-10 10])
hold on;
plot(t2, ysolution2(:,1), 'r');

ylabel('x');
xlabel('t');
%axis('tight', 'equal')
subplot(1,2,2);
plot(t, ysolution(:,2), 'b');
hold on;
plot(t2, ysolution2(:,2), 'r');
legend('x_0,y_0 = (0, -3)', 'x_0,y_0 = (25, 3)')
ylabel('y')
xlabel('t');

hold on;
%% Part 2 - Homoclinic Bifurcation
mu_values = [-0.5 -.2 0 .2 .5 1];
figure();
for i=1:length(mu_values)
    subplot(3,2,i)
    mu = mu_values(i);
    
    dx_dt = @(x,y) mu*x+y-x.^2;
    dy_dt = @(x,y) -x+mu*y+2*x.^2;
    x_range = -2:0.02:2;
    nullcline_x = x_range.^2 - mu*x_range;
    nullcline_y = (1/mu)*(x_range-2*x_range.^2);
    [x, y] = meshgrid(-2:0.02:2, -2:0.02:2);
    
    dx = dx_dt(x, y);
    dy = dy_dt(x, y);
    
    streamslice(x, y, dx, dy);
    hold on;
    plot(x_range, nullcline_x, 'k--');
    hold on;
    plot(x_range, nullcline_y, 'k--');
    xlim([-2 2])
    ylim([-2 2])
    xlabel('x'), ylabel('y')
    title(strcat('\mu=', num2str(mu)))
    %axis('tight', 'equal')
    hold on;
end
%% %% Part 2 - Homoclinic Bifurcation
mu_values = [.015 .02 .025 .05 .065 .07];
figure();
for i=1:length(mu_values)
    subplot(3,2,i)
    mu = mu_values(i);
    
    dx_dt = @(x,y) mu*x+y-x.^2;
    dy_dt = @(x,y) -x+mu*y+2*x.^2;
    x_range = -2:0.02:2;
    nullcline_x = x_range.^2 - mu*x_range;
    nullcline_y = (1/mu)*(x_range-2*x_range.^2);
    [x, y] = meshgrid(-2:0.02:2, -2:0.02:2);
    
    dx = dx_dt(x, y);
    dy = dy_dt(x, y);
    
    streamslice(x, y, dx, dy);
    hold on;
    plot(x_range, nullcline_x, 'k--');
    hold on;
    plot(x_range, nullcline_y, 'k--');
    xlim([-0.75 0.75])
    ylim([-0.75 0.75])
    xlabel('x'), ylabel('y')
    title(strcat('\mu=', num2str(mu)))
    %axis('tight', 'equal')
    hold on;
end
sgtitle('Homoclinic Bifurcation: Phase Portraits')
%% 
options = odeset('RelTol',1e-8);
mu_values = [.015 .02 .025 .05 .065 .07];
figure();

for i=1:length(mu_values)
    subplot(2,1,1)
    mu = mu_values(i);
    
    [t,ysolution] = ode45(@(t,y) [mu*y(1)+y(2)-y(1).^2; -1*y(1)+mu*y(2) + 2*y(1).^2],[0 45],[-0.05 -0.05],options);
    plot(t, ysolution(:,1));
    ylabel('x')
    ylim([-0.5 0.5])
    hold on;
    subplot(2,1,2)
    plot(t, ysolution(:,2));
    ylim([-0.5 0.5])
    xlabel('t'), ylabel('y')
    %axis('tight', 'equal')
    hold on;
end
sgtitle('Homoclinic Bifurcation: Time Traces')
subplot(2,1,1);
legend(strcat('\mu=',string(num2cell(mu_values))))


%% 
options = odeset('RelTol',1e-8);
mu_values = [.065 .067 .07];
initial_values = -.1:.05:.1;
figure();

for i=1:length(mu_values)
    subplot(3,2,i)
    mu = mu_values(i);
    for j=1:length(initial_values)
        [t,ysolution] = ode45(@(t,y) [mu*y(1)+y(2)-y(1).^2; -1*y(1)+mu*y(2) + 2*y(1).^2],[0 50],[initial_values(j) initial_values(j)],options);
        plot( ysolution(:,1), ysolution(:,2));
        ylabel('x')
        ylim([-0.75 0.75])
        ylim([-0.75 0.75])
        xlabel('x'), ylabel('y')
        %axis('tight', 'equal')
        hold on;
        title(strcat('\mu=', num2str(mu)))
    end
end
sgtitle('Homoclinic Bifurcation')
