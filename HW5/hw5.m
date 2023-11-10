%% 
mu = 4;
figure();
[x, y] = meshgrid(-10:0.02:10, -10:0.02:10);
    
dx_dt = @(x,y) mu*x - x.^3;
dy_dt = @(x,y) -y;
dx = dx_dt(x, y);
dy = dy_dt(x, y);

streamslice(x, y, dx, dy);
%% 
% Problem 3
mu_values = [-2 -.2 0 .1 .2 ];
figure();
for i=1:length(mu_values)
    subplot(3,2,i)
    mu = mu_values(i);
    
    x_range = -2:0.02:2;
    % nullcline_x = x_range.^2 - mu*x_range;
    % nullcline_y = (1/mu)*(x_range-2*x_range.^2);
    [x, y] = meshgrid(-2:0.02:2, -2:0.02:2);
    
    dx_dt = @(x,y) -y + mu*x +x.*y.^2;
    dy_dt = @(x,y) x + mu*y - x.^2;
    dx = dx_dt(x, y);
    dy = dy_dt(x, y);

    streamslice(x, y, dx, dy);
    hold on;
    % plot(x_range, nullcline_x, 'k--');
    % hold on;
    % plot(x_range, nullcline_y, 'k--');
    xlim([-2 2])
    ylim([-2 2])
    xlabel('x'), ylabel('y')
    title(strcat('\mu=', num2str(mu)))
    %axis('tight', 'equal')
    hold on;
end

hold on;

%% 
options = odeset('RelTol',1e-8);
% mu_values = [-2 -.2 0 .1 .2 ];
% mu_values = [-2 -.01 0];
mu_values = [-.1]

% figure();
for i=1:length(mu_values)
    figure();
    
    % subplot(2,1,1)
    mu = mu_values(i);
    
    [t,ysolution] = ode45(@(t,y) [-y(2) + mu*y(1) + y(1).*y(2).^2; ...
        y(1)+mu*y(2)-y(1).^2],[0 3000],[-.1 -.1],options);
    % plot(t, ysolution(:,1));
    % plot(ysolution(:,1), ysolution(:,2))
    title(strcat('\mu=', num2str(mu)))
    hold on;
    ylabel('x')
    % ylim([-10 10])
    hold on;
    % title(strcat('\mu=', num2str(mu)))
    % subplot(2,1,2)
    % plot(t, ysolution(:,2));
    % ylim([-10 10])

    
    xlabel('x'), ylabel('y')
    [t,ysolution] = ode45(@(t,y) [-y(2) + mu*y(1) + y(1).*y(2).^2; ...
        y(1)+mu*y(2)-y(1).^2],[0 3000],[.0001 .0001],options);
    % subplot(2,1,2)
    plot(ysolution(:,1), ysolution(:,2))
    % legend('small', 'very small')
    %axis('tight', 'equal')
    hold on;
    

end
% sgtitle('Hopf Bifurcation: Phase Portraits')
% subplot(2,1,1);
% legend(strcat('\mu=',string(num2cell(mu_values))))

%% 
options = odeset('RelTol',1e-8);
mu_values = [-2 -.2 0 .1 .2 ];
figure();

for i=1:length(mu_values)
    initial_values = -1:.1:1;
    subplot(3,2,i)
    mu = mu_values(i);
    for j=1:length(initial_values)
        [t,ysolution] = ode45(@(t,y) [-y(2) + mu*y(1) + y(1).*y(2).^2; ...
            y(1)+mu*y(2)-y(1).^2],[0 45],[initial_values(j) initial_values(j)],options);
        plot(ysolution(:,1), ysolution(:,2));
        title(strcat('\mu=', num2str(mu)))
        hold on;
    
        ylabel('x')
        ylim([-1 1])
        hold on;
        dx_dt = @(x,y) -y + mu*x +x.*y.^2;
        dy_dt = @(x,y) x + mu*y - x.^2;
        dx = dx_dt(x, y);
        dy = dy_dt(x, y);

    streamslice(x, y, dx, dy, '--');
        % subplot(2,1,2)
        % plot(t, ysolution(:,2));
        % ylim([-.5 .5])
        % xlabel('t'), ylabel('y')
        % %axis('tight', 'equal')
        % hold on;
    end
end
% sgtitle('Hopf Bifurcation: Time Traces')
% subplot(2,1,1);
% legend(strcat('\mu=',string(num2cell(mu_values))))

%% 
% mu_values = [-1 -.2 0 .1 .2 ];
mu=2;
figure();
for i=1:length(mu_values)
    % mu = mu_values(i);
    
    [t,ysolution] = ode45(@(t,y) [y(2)-2*(y(1)); mu + y(1).^2 - y(2)],[0 100],[1 1],options);
    plot(t, ysolution(:,1));
    ylabel('x')
    hold on;
    title(strcat('\mu=', num2str(mu)))
    subplot(2,1,2)
    ylim([-20 20])
    plot(t, ysolution(:,2));
    %axis('tight', 'equal')
    hold on;
end

hold on;
