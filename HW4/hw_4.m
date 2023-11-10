%% 
params = {};
params.g_in = 0;
params.k1 = 0.5;
params.k2 = 0.1;
params.k3 = 1.0;
params.k4 = 0.1;
params.tau = 15;
params.n = 4;
%% test a bunch of initial conditions
t_end = 500;

opts = ddeset('RelTol',1e-8,'AbsTol',1e-8);



initial_x = 1:1:5;
initial_y =  1:1:5;
figure();
for i=1:length(initial_x)
    for j=1:length(initial_y)
        sol = dde23(@(t,x,Z) gi_eq(t,x,Z,params),params.tau,@(t) ...
                [initial_x(i) initial_y(j)],[0, t_end],opts);

        subplot(2,1,1);
        plot(sol.x,sol.y(1,:),'--','LineWidth',1);
        ylabel('Glucose (G)')
        hold on;
        subplot(2,1,2);
        plot(sol.x,sol.y(2,:),'--','LineWidth',1);
        ylabel('Insulin (I)')
        hold on;
    end
end
sgtitle('Glucose/Insulin Dynamics')
xlabel('Time (min)')
%% Just one trace
params.tau = 15;
params.n = 4;
figure();

sol = dde23(@(t,x,Z) gi_eq(t,x,Z,params),params.tau,@(t) ...
        [2 2],[0, t_end],opts);

subplot(2,1,1);
plot(sol.x,sol.y(1,:),'--','LineWidth',1);
hold on;
[peaks, indices] = findpeaks(sol.y(1,:));
plot(sol.x(indices), peaks, 'o');

ylabel('Glucose (G)')

subplot(2,1,2);
plot(sol.x,sol.y(2,:),'--','LineWidth',1);

ylabel('Insulin (I)')
hold on;
[peaks, indices] = findpeaks(sol.y(2,:));
plot(sol.x(indices), peaks, 'o');

sgtitle('Glucose/Insulin Dynamics')
xlabel('Time (min)')
%% 
tau_values = 1:1:4;
n_values = 1:1:8;
for t=1:length(tau_values)
    for n=1:length(n_values)
        params.tau = tau_values(t);
        params.n = n_values(n);
        figure();
        
        sol = dde23(@(t,x,Z) gi_eq(t,x,Z,params),params.tau,@(t) ...
                [2 2],[0, t_end*2],opts);
        
        
        subplot(2,1,1);
        plot(sol.x,sol.y(1,:),'--','LineWidth',1);
        hold on;
        
        ylabel('Glucose (G)')
        
        subplot(2,1,2);
        plot(sol.x,sol.y(2,:),'--','LineWidth',1);
        
        ylabel('Insulin (I)')
        hold on;
        
        
        xlabel('Time (min)')

        title(strcat('tau=' ,num2str(params.tau), '; n=', num2str(params.n)))


    end
end

%%  ntau plot
figure()

hold on;

scatter([2 3 4 5 6 7 8], 15, 100,'r', 'filled');
scatter(1, 15, 100,'b', 'filled');

scatter([1 2 3 4 5 6 7 8], 20, 100,'r', 'filled');
scatter([2 3 4 5 6 7 8], 10, 100, 'r', 'filled');
scatter(1, 10, 100, 'b', 'filled');
scatter([3 4 5 6 7 8], 5, 100, 'r', 'filled');
scatter([1 2], 5, 100, 'b', 'filled')
scatter([3 4 5 6 7 8], 4, 100, 'r', 'filled');
scatter([1 2], 4, 100, 'b', 'filled');

scatter([4 5 6 7 8], 3, 100, 'r', 'filled');
scatter([1 2 3], 3, 100, 'b', 'filled');

scatter([6 7 8], 2, 100, 'r', 'filled');
scatter([1 2 3 4 5], 2, 100, 'b', 'filled');
scatter([1 2 3 4 5 6 7 8], 1, 100, 'b', 'filled');
xlabel('n')
ylabel('Tau')
xlim([0 11])
ylim([0 21])
text(9, 20, {'{\color{red} o } Stable Oscillations', '{\color{blue} o } Non-Oscillatory'});
%% 
p = 0:.1:10;
s = 0:.1:10;
figure();
plot(p, 1./(.9*p.*p));
hold on;

plot([0 0], [0 5]);
plot(p, 1./(.9*p));
plot(1,1/.9, 'o')
xlim([-1 6])
ylim([0, 5])
legend('(1)', '(2)', '(3)', 'Fixed Point')
xlabel('P')
ylabel('S')
title('Nullclines')

[t,ysolution] = ode45(@(t,y) [c.*y(2).*y(1).*y(1) - k*y(1); ...
            v0-c.*y(2).*y(1).*y(1)],[0 t_end],[1, 12],options);
plot(ysolution(:,1),ysolution(:,2),'--','LineWidth',0.5);


%% 
% [pgrid,sgrid] = meshgrid(0:0.1:1.8,0:0.1:1.8);
% c = 0.9;
% k=1;
% v0=1;
% pdot = c*sgrid.*pgrid.*pgrid - k*pgrid;
% sdot = v0 - c*sgrid.*pgrid.*pgrid;
% figure();
% subplot(3,1,1);
% plot(1,10/9, 'o')
% hold on;
% quiver(pgrid, sgrid, pdot, sdot);
hold on;
p_0 = 0.01:0.3:.8;
s_0 = 0.01:0.3:.8;
t_end = 100;

options = odeset('RelTol',1e-8);
for p=1:length(p_0)
    for s=1:length(s_0)
        [t,ysolution] = ode45(@(t,y) [c.*y(2).*y(1).*y(1) - k*y(1); ...
            v0-c.*y(2).*y(1).*y(1)],[0 t_end],[p_0(p), s_0(s)],options);
        % figure();
        subplot(3,1,1);
        
        plot(ysolution(:,1),ysolution(:,2),'--','LineWidth',0.5);
        hold on;
        xlabel('P')
        ylabel('S')
        subplot(3,1,2);
        hold on;
        plot(t, ysolution(:,1));
        xlabel('t')
        ylabel('P')
        subplot(3,1,3);
        hold on;
        plot(t, ysolution(:,2));
        xlabel('t')
        ylabel('S')

    end
end
subplot(3,1,1);


%% 

function dydt = diff_eq(t,y,Z)
% Differential equations function 
g_in = 0;
k1 = 0.5;
k2 = 0.1;
k3 = 1.0;
k4 = 0.1;
tau = 15;
n = 4;

ylag = Z(:,1);
dydt = [ g_in + k3/(1+y(2)^2) - k4*y(1) - y(1).*y(2);
         (k1*ylag(1).^n)/(1+ylag(1).^n) - k2*(y(2))];

end

function dydt = gi_eq(t,y,Z,params)

% Differential equations function     
ylag = Z(:,1);
dydt = [ params.g_in + params.k3/(1+y(2)^2) - params.k4*y(1) - y(1).*y(2);
         (params.k1*ylag(1).^params.n)/(1+ylag(1).^params.n) - params.k2*(y(2))];

end