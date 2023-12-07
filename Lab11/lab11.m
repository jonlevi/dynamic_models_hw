global L;
global t_end;
global D;
global x_on;
global V_on;
global gbar_K;
global gbar_Ca;

L = 30; %length of neuron (in cm)
t_end = 300; %time to simulate over (in ms)
D = 0.0375; %diffusion coeffcient cm^2/ms
x_on = 1; %cm, length over which V is elevated initially
V_on = 20; %mV; value V is elevated to initially
gbar_K=8;
gbar_Ca=4.4;

% gca_values = 4.5:.5:6;


% starting_voltages = [-30 -28 -27  -25];
% D_values = .01:.01:.1 ;
speeds = zeros(1, length(gca_values));
for s=1:1

    % 
    % V_on = starting_voltages(s);
    % D = D_values(s);
    % gbar_K = gk_values(s);
    % gbar_Ca = gca_values(s);

    %set up solution mesh
    m = 0;
    Nx = 200;
    Nt = 300;
    x = linspace(0,L,Nx);
    t = linspace(0,t_end,Nt);
    
    %solve the partial differential equations
    sol = pdepe(m,@Morris_Lecar_cable_eqn,@Morris_Lecar_cable_initial,@Morris_Lecar_cable_bc,x,t);
    
    V = sol(:,:,1);
    w = sol(:,:,2);
    
    % % plot results
    figure(1)
    % subplot(2,2,s);
    surf(x,t,V, 'EdgeColor','none')    
    xlabel('Position (cm)')
    ylabel('Time (ms)')
    zlabel('Membrane Voltage (mV)')
    shading interp
    % title(strcat('gca= ', num2str(gbar_Ca)))
    set(gca,'fontsize',14)
    % view(2)
    colorbar;
    clim([-60, 35]);
    hold on;
    [~, index] = max(V(200,:));
    % plot([x(index) x(index)] ,[0 200], 'w--')
    speeds(s) = x(index)/200;

    
    % 
    % %plot spatial profiles at a few different time points
    figure(2)
    subplot(2,1,1);
    set(gca,'fontsize',14)
    hold on
    plot(x, V(1,:),'k', 'linewidth', 2);
    plot(x, V(50,:),'k--', 'linewidth', 2);
    plot(x, V(150,:),'k:', 'linewidth', 2);
    %axis([0 15 -80 40])
    xlabel('Position (cm)')
    ylabel('Membrane voltage (mV)')
    title(strcat('Starting Voltage: ', num2str(V_on)))
    legend('t=0 ms', 't=50 ms', 't=150 ms')

    subplot(2,1,2);
    set(gca,'fontsize',14)
    hold on
    plot(x, V(200,:),'k', 'linewidth', 2);
    plot(x, V(250,:),'k--', 'linewidth', 2);
    plot(x, V(300,:),'k:', 'linewidth', 2);
    %axis([0 15 -80 40])
    xlabel('Position (cm)')
    ylabel('Membrane voltage (mV)')
    title(strcat('Starting Voltage: ', num2str(V_on)))
    legend('t=200 ms', 't=250 ms', 't=300 ms')

    % figure(3);
    % plot(t, V(:,1));
    % xlabel('time')
    % ylabel('V(t) [at x=1]')
    % hold on;
end
% figure(3);
% legend(strcat('V0=',string(num2cell(starting_voltages))));
% figure();
% plot(gca_values, speeds);
% xlabel('gCa bar')
% ylabel('Speed (cm/ms)');
%% 


%set equations
function [c,f,s] = Morris_Lecar_cable_eqn(x,t,u,DuDx)

global D;
global gbar_Ca;
global gbar_K;

Cap=20;
V_Nernst_K=-84;
V_Nernst_Ca=120;
gbar_leak=0.5;
V_Nernst_leak=-60;
v1=-1.2;
v2=18;
v3=2;
v4=30;

m_inf = @(x) ( 0.5*(1+tanh((x-v1)/v2)) );
w_inf = @(x) ( 0.5*(1+tanh((x-v3)/v4)) );
tau_w=0.8/.04;

c = [1;1];
f = [D;0].*DuDx;
s = [(1/Cap)*(-gbar_Ca*m_inf(u(1))*(u(1)-V_Nernst_Ca)-gbar_K*u(2)*(u(1)-V_Nernst_K)...
        -gbar_leak*(u(1)-V_Nernst_leak)); 
    (w_inf(u(1))-u(2))/tau_w];

end


%set initial condition
function value = Morris_Lecar_cable_initial(x)

global x_on;
global V_on;
global L;

if (x<=x_on) || (x >= L-x_on)
    value = [V_on; 0.004];
else
    value = [-61.9; 0.004];
end

end


%set boundary conditions
% of the form p + q (D dv/dx)) = 0
function [pl,ql,pr,qr] = Morris_Lecar_cable_bc(xl,ul,xr,ur,t)

% Left side
pl = 0;
ql = 1;

% Right side
a = 1;
pr = 
qr = [1;1];

end