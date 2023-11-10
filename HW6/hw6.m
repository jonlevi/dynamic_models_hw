% Lorenz Model

options = odeset('RelTol',1e-8);
t_end = 50;
model = @(t,y) [10*(y(2)-y(1)); (y(1).*(28-y(3)))-y(2); y(1).*y(2)-((8/3)*y(3))];
x_range = -20:20:40;
y_range = -25:15:25;
z_range = 0:20:50;
for i=1:length(x_range)
    for j=1:length(y_range)
        for q=1:length(z_range)
            y0 = [x_range(i) y_range(j) z_range(q)];
            
            [t,y] = ode45(model,[0 t_end],y0,options);
            figure(1)
            subplot(3,1,1);
            plot(t,y(:,1),'-','LineWidth',1);
            hold on
            ylabel('x');
            set(gca,'Fontsize',14,'LineWidth',1)
            subplot(3,1,2);
            plot(t,y(:,2),'-','LineWidth',1)
            hold on
            set(gca,'Fontsize',14,'LineWidth',1);
            ylabel('y');
            subplot(3,1,3);
            plot(t,y(:,3),'-','LineWidth',1);
            hold on;
            ylabel('z');
            xlabel('Time (t.u.)'), 
            
            figure(2);
            plot3(y(:,1), y(:,2), y(:,3));
            xlabel('x')
            ylabel('y')
            zlabel('z')
            hold on;
        end
    end
end
%% 

% simulate starting on attractor by taking the end of a diff simulation
[t,y] = ode45(model,[0 t_end],[-20 -10 20],options);
figure(3);
plot3(y(:,1), y(:,2), y(:,3));
hold on;
N = length(t);
plot3(y(N,1), y(N,2), y(N,3),  '.r', 'MarkerSize', 30);
legend('',"(" + num2str(y(N,1)) + ", " +  num2str(y(N,2)) + ", " + num2str(y(N,3)) + ")")
xlabel('x')
ylabel('y')
zlabel('z')

% using end as new starting point
y0 = [y(N,1), y(N,2), y(N,3)];
[t,y2] = ode45(model,[0 t_end],y0,options);
figure(1)
subplot(3,1,1);
plot(t,y2(:,1),'-','LineWidth',1);
hold on
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
subplot(3,1,2);
plot(t,y2(:,2),'-','LineWidth',1)
hold on
set(gca,'Fontsize',14,'LineWidth',1);
ylabel('y');
subplot(3,1,3);
plot(t,y2(:,3),'-','LineWidth',1);
hold on;
ylabel('z');
xlabel('Time (t.u.)'), 

figure(2);
plot3(y2(:,1), y2(:,2), y2(:,3));
hold on;
plot3(y(N,1), y(N,2), y(N,3),  '.r', 'MarkerSize', 30);
legend('',"(" + num2str(y(N,1)) + ", " +  num2str(y(N,2)) + ", " + num2str(y(N,3)) + ")")
xlabel('x')
ylabel('y')
zlabel('z')
hold on;

%% 
% initial condition inside attractor
t_end = 15;
y0 = [y(N,1), y(N,2), y(N,3)];
[t2,y2] = ode45(model,[0 t_end],y0,options);
figure();
delta_values = [.001, .01, .05, .1];
for d=1:length(delta_values)
    subplot(2,2,d);
    plot(t2,y2(:,1),'-','LineWidth',1);
    hold on;
    delta = delta_values(d);
    y0_delta = [y(N,1)+delta, y(N,2)+delta, y(N,3)+delta];
    [t3,y3] = ode45(model,[0 t_end],y0_delta,options);
    plot(t3,y3(:,1),'-','LineWidth',1);
    title(strcat('\Delta=',num2str(delta)))
   
end
subplot(2,2,2);
legend('Y0', 'Y0 + \Delta');
