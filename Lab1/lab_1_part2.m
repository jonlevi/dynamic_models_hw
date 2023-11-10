% Lactose Switch Model (Lab 1 Part 2)

% initial condition
x0 = 5;

%parameters
t_end = 50;
k = 0.3;
a = 0.05;

dxdt = @(t,x) (a+x.^2)/(1+x.^2) - (k*x);

options = odeset('RelTol',1e-8);
[t,x] = ode45(dxdt,[0 t_end],x0,options);

%plot numerical solution
hold on
plot(t,x,'b.-')
xlabel('t'), ylabel('x')
title('Lactose switch model')
hold off

%% try changing initial conditions with same parameters
initial_conditions = 0:.5:5;
k = 0.3;
dxdt = @(t,x) (a+x.^2)/(1+x.^2) - (k*x);
figure()

for i = 1:length(initial_conditions)
    [t,x] = ode45(dxdt,[0 t_end],initial_conditions(i),options);
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
end
title('Lactose switch model, k=0.3')
legend(strcat('X0=',string(num2cell(initial_conditions))))
%%

% try changing initial conditions with k=0.6
hold off
figure()
k = 0.6;
dxdt = @(t,x) (a+x.^2)/(1+x.^2) - (k*x);
for i = 1:length(initial_conditions)
    [t,x] = ode45(dxdt,[0 t_end],initial_conditions(i),options);
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
end
legend(strcat('X0=',string(num2cell(initial_conditions))))
title('Lactose switch model, k=0.6')
hold off


%%
%Bistability anaylsis
a = 0.05;
figure();
hold on;
x = linspace(0,5, 500);
g = (a+x.^2)./(1+x.^2);
plot(x, g, 'r-');

k_values=linspace(.1,1,10);
for i=1:length(k_values)
    h = k_values(i)*x;
    plot(x, h, '--');
    title('Components of the Lac Switch')
end
legend(["Import" strcat('Metabolism k=',string(num2cell(k_values)))])
hold off
%%

% You can see that the curve with k=0.5 has three intersections with g
% and is thus bistable
%%
k = 0.5;
dxdt = @(t,x) (a+x.^2)/(1+x.^2) - (k*x);
figure()
hold on;
initial_conditions = linspace(0,2,15);
for i = 1:length(initial_conditions)
    [t,x] = ode45(dxdt,[0 t_end],initial_conditions(i),options);
    %plot numerical solution
    hold on
    plot(t,x)
    xlabel('t'), ylabel('x')
end
legend(strcat('X0=',string(num2cell(initial_conditions))))
title('Lactose switch model, k=0.5')




%%
k_values = .42:.01:.54;
figure()
hold on;
for j=1:length(k_values)
    dxdt = @(t,x) (a+x.^2)/(1+x.^2) - (k_values(j)*x);
    subplot(7,2,j)
    initial_conditions = linspace(0,2,15);
    for i = 1:length(initial_conditions)
        [t,x] = ode45(dxdt,[0 t_end*15],initial_conditions(i),options);
        %plot numerical solution
        hold on
        plot(t,x)
        xlabel('t'), ylabel('x')
    end
    title(strcat('k=', num2str(k_values(j))))
    
end


