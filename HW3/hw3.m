% HW 3
[xgrid,ygrid] = meshgrid(-10:1:20,-5:1:10);
xdot = ygrid + ygrid.*ygrid;
ydot = -(1/2)*xgrid + (1/5)*ygrid - xgrid.*ygrid + (6/5)*ygrid.*ygrid;
figure()
quiver(xgrid, ygrid, xdot, ydot)
hold on;
t_end = 5;
options = odeset('RelTol',1e-8);

% initial conditions
x_0 = -2:.5:2;
y_0 = -5:.5:2;
for i=1:length(x_0)
    for j=1:length(y_0)
        [t,ysolution] = ode45(@(t,x) [x(2) + x(2).*x(2); ...
            -(1/2)*x(1) + (1/5)*x(2)- x(1).*x(2) + (6/5)*x(2).*x(2)], ...
            [0 t_end],[x_0(i), y_0(j)],options);

        %plot solution
        plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
        hold on;
        set(gca,'Fontsize',14,'LineWidth',1)

    end
end
xlabel('X')
ylabel('Y')
title('Phase Portrait')

%% 

[xgrid,ygrid] = meshgrid(-3:.1:3,-3:.1:3);
xdot = xgrid - ygrid;
ydot = (xgrid.*xgrid) -4 ;
figure()
quiver(xgrid, ygrid, xdot, ydot)
hold on;
plot(-3:.1:3, -3:.1:3, 'k--')
plot([-2 -2], [-3, 3], 'k--')
plot([2 2], [-3, 3], 'k--')
plot(-2,-2, 'ko');
plot(2,2, 'ko');
t_end = 1.2;
options = odeset('RelTol',1e-8);

% initial conditions
% x_0 = -.05:.01:.05;
% y_0 = -.05:.01:.05;
% x_0 = -0.055:.01:.02;
x_0 = -3:.5:3;
y_0 = -3:.5:3;
for i=1:length(x_0)
    for j=1:length(y_0)
    
        [t,ysolution] = ode45(@(t,x) [x(1)-x(2); x(1).*x(1)-4], ...
            [0 t_end],[x_0(i), y_0(j)],options);
        
        %plot solution
        plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
        hold on;

        xlim([-3.2 3.2]);
        ylim([-3.2 3.2]);
    end
end
% plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
% [t,ysolution] = ode45(@(t,x) [x(1)-x(2); x(1).*x(1)-4], ...
%     [0 t_end],[0, -1],options);
% plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
% [t,ysolution] = ode45(@(t,x) [x(1)-x(2); x(1).*x(1)-4], ...
%     [0 t_end],[-2, -1.8],options);
% plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
% [t,ysolution] = ode45(@(t,x) [x(1)-x(2); x(1).*x(1)-4], ...
%     [0 t_end],[-1.9, .2],options);
% plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
% set(gca,'Fontsize',14,'LineWidth',1)


xlabel('X')
ylabel('Y')
title('Phase Portrait')
%% 
[xgrid,ygrid] = meshgrid(-2:.1:2,-1.2:.1:1.2);
t_end = 10;
xdot = ygrid + xgrid - xgrid.^3;
ydot = -ygrid;
samplex = -2:.1:2;
figure()
quiver(xgrid, ygrid, xdot, ydot)
hold on;
plot([-2 2], [0 0], 'k--')
plot(samplex, samplex.^3-samplex, 'k--')
plot(0,0, 'ko');
plot(-1,0, 'ko');
plot(1,0, 'ko');
ylim([-1.2 1.2])
x_0 = -2:.5:2;
y_0 = -1:.5:1;
for i=1:length(x_0)
    for j=1:length(y_0)
    
        [t,ysolution] = ode45(@(t,x) [x(2)+x(1)-x(1).^3; -x(2)], ...
            [0 t_end],[x_0(i), y_0(j)],options);
        
        %plot solution
        plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
        hold on;
    end
end
%% Bonus problem
x = 
xdot = r + x.^2;
