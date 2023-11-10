%% 
% initial conditions
x0_values = -5:1:3;

%parameters
t_end = .2;

for i=1:length(x0_values)

    options = odeset('RelTol',1e-8);
    [t,x] = ode45(@(t,x) 4*x.^2-16,[0 t_end],x0_values(i),options);
    
    
    dt = t_end/(length(t)-1);
    t_even = 0:dt:t_end;
    
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
    title('Solving dx/dt = 4x^2-16')
    ylim([-5 5])
end
legend(strcat('X0=',string(num2cell(x0_values))))

%% 
figure()
x0_values = -1.5:.5:3;

%parameters
t_end = 3;

for i=1:length(x0_values)

    options = odeset('RelTol',1e-8);
    [t,x] = ode45(@(t,x) 1-x.^14,[0 t_end],x0_values(i),options);
    
    
    dt = t_end/(length(t)-1);
    t_even = 0:dt:t_end;
    
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
    title('Solving dx/dt = 1-x^{14}')
    xlim([-.1 3])

end
legend(strcat('X0=',string(num2cell(x0_values))))


%% 
figure()
x0_values = -3:.5:3;

%parameters
t_end = 5;

for i=1:length(x0_values)

    options = odeset('RelTol',1e-8);
    [t,x] = ode45(@(t,x) x-x.^3,[0 t_end],x0_values(i),options);
    
    
    dt = t_end/(length(t)-1);
    t_even = 0:dt:t_end;
    
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
    title('Solving dx/dt = x-x^{3}')
end
legend(strcat('X0=',string(num2cell(x0_values))))


%% 
figure()
x0_values = 2*pi:.5:3*pi;

%parameters
t_end = 50000;

for i=1:length(x0_values)

    options = odeset('RelTol',1e-8);
    [t,x] = ode45(@(t,x) exp(-x).*sin(x),[0 t_end],x0_values(i),options);
    
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
    title('Solving dx/dt = e^{-x}sinx')
    xlim([-1 t_end])
end
legend(strcat('X0=',string(num2cell(x0_values))))
%% 
figure()
x0_values = -2:.5:2;

%parameters
t_end = 50;

for i=1:length(x0_values)

    options = odeset('RelTol',1e-8);
    [t,x] = ode45(@(t,x) 1+.5*cos(x),[0 t_end],x0_values(i),options);
    
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
    title('Solving dx/dt = 1+ 1/2 cos(x)')
    xlim([-1 t_end])
end
legend(strcat('X0=',string(num2cell(x0_values))))