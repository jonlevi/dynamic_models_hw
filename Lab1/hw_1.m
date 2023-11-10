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
