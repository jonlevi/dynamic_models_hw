%% quiver plots part 1
a_values = [-1.5 -1 -0.5 0 1];
figure();
hold on;
for i=1:length(a_values)
    a = a_values(i);
    [x,y] = meshgrid(-2:0.2:2,-2:0.2:2);
    xdot = a*x;
    ydot = -y;
    subplot(3,2,i)
    quiver(x,y,xdot,ydot)
    title('A='+ string(a))
end

%% simulate solutions from first part
a_values = [-1.5 -1 -0.5 0 1];
t_end = 5;
options = odeset('RelTol',1e-8);

for i=1:length(a_values)
    figure();
    hold on;
    a = a_values(i);
    % exponential growth makes it hard to visualize
    if a > 0
        t_end = 1;
    end
    initial_x = -2:0.2:2;
    initial_y = -2:0.2:2;
    for x_index=1:length(initial_x)
        for y_index=1:length(initial_y)
            x_0 = initial_x(x_index);
            y_0 = initial_y(y_index);


        
            [t,x] = ode45(@(t,x) [a*x(1);-x(2)],[0 t_end],[x_0; y_0],options);
        
            [xgrid,ygrid] = meshgrid(-2:0.2:2,-2:0.2:2);
            subplot(3,1,1)
            plot(t,x(:,1))
            title('x(t)')
            subplot(3,1,2)
            plot(t,x(:,2))
            title('y(t)')
            subplot(3,1,3)
            plot(x(:,1),x(:,2), 'b')
            hold on;
            quiver(xgrid,ygrid,a*xgrid,-1*ygrid, 'r')
            title('Vector Plots (a=' + string(a) + ')')
        end
    end
    print('-dpdf', 'Lab2/figure' + string(i+1) + '.pdf', '-fillpage')
end
%% 

%% quiver plot part 2
figure();
hold on;
[x,y] = meshgrid(-2:0.2:2,-2:0.2:2);
xdot = y;
ydot = -x+y;
quiver(x,y,xdot,ydot)
title('Quiver Plot for Love Affair Model')
print('-dpdf', 'Lab2/figure7.pdf', '-fillpage')

%% 
initial_x = -1:0.2:1;
initial_y = -1:0.2:1;
t_end = 5;
figure();
hold on;
for x_index=1:length(initial_x)
    for y_index=1:length(initial_y)
        x_0 = initial_x(x_index);
        y_0 = initial_y(y_index);
        [xgrid,ygrid] = meshgrid(-10:1.5:10,-10:1.5:10);
        xdot = ygrid;
        ydot = -xgrid+ygrid;
        [t,x] = ode45(@(t,x) [x(2);-x(1)+x(2)],[0 t_end],[x_0; y_0],options);
        subplot(3,1,1)
        plot(t,x(:,1))
        title('x(t)')
        subplot(3,1,2)
        plot(t,x(:,2))
        title('y(t)')
        subplot(3,1,3)
        plot(x(:,1),x(:,2), 'b')
        xlim([-15 15])
        hold on;
        quiver(xgrid,ygrid,xdot,ydot, 3, 'r')
        title('Love Affair Model')
    end
end

