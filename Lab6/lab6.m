
% initial conditions
y0 = [-30; -30];

%parameters
t_end = 100;

options = odeset('RelTol',1e-8,'Maxstep',0.1);
[t,y] = ode15s(@FHNeqns,[0 t_end],y0,options);

%plot numerical solution
figure(1)
subplot(211), plot(t,y(:,1),'b-','LineWidth',2)
hold on
ylabel('x')
set(gca,'Fontsize',14,'LineWidth',1)
subplot(212), plot(t,y(:,2),'b-','LineWidth',2)
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')

%plot numerical solution
figure(2)
plot(y(:,1),y(:,2),'b-','LineWidth',2)
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('x'), ylabel('y')

%% 
t_end = 100;
options = odeset('RelTol',1e-8,'Maxstep',0.1);
initial_x = -1:.5:1;
initial_y = -1:.5:1;
for i=1:length(initial_x)
    for j=1:length(initial_y)
        y0 = [initial_x(i); initial_y(j)];
        [t,y] = ode15s(@FHNeqns,[0 t_end],y0,options);
        %plot numerical solution
        figure(1)
        subplot(211), plot(t,y(:,1),'--','LineWidth',.5)
        hold on
        ylabel('x')
        set(gca,'Fontsize',14,'LineWidth',1)
        subplot(212), plot(t,y(:,2),'--','LineWidth',.5)
        hold on
        set(gca,'Fontsize',14,'LineWidth',1)
        xlabel('Time (t.u.)'), ylabel('y')
        
        %plot numerical solution
        figure(2)
        plot(y(:,1),y(:,2),'--','LineWidth',.5)
        hold on
        set(gca,'Fontsize',14,'LineWidth',1)
        xlabel('x'), ylabel('y')

    end
end

%% 


% initial conditions
y0 = [-1.2; -0.62];

%parameters
t_end = 100;

options = odeset('RelTol',1e-8,'Maxstep',0.1);
[t,y] = ode15s(@(t,y) FHNExcitability(t,y,.459),[0 t_end],y0,options);
[t2,y2] = ode15s(@(t,y) FHNExcitability(t,y,0.46),[0 t_end],y0,options);
%plot numerical solution
figure(1)
subplot(2,1,1);
plot(t,y(:,1),'b-');
hold on;
plot(t2,y2(:,1), 'r-');
legend('Amp=0.459', 'Amp=0.46')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
subplot(2,1,2);
plot(t,y(:,2),'b-');
hold on;
plot(t2,y2(:,2), 'r-');
hold on;

set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')
legend('Amp=0.459', 'Amp=0.46')
title('Excitability')

%plot numerical solution
figure(2)
plot(y(:,1),y(:,2),'b-','LineWidth',2)
hold on;
plot(y2(:,1),y2(:,2),'r-','LineWidth',2)

x_range = -2.5:.1:2;
nullcline_1 = (1/0.8)*(x_range+0.7);
nullcline_2 = (-1/3)*x_range.^3+x_range;
plot(x_range, nullcline_1, 'g--')
plot(x_range, nullcline_2, 'k--')

legend('Amp=1.0', 'Amp=0.3', 'Ydot nullcline', 'Xdot nullcline')
title('Excitability')
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('x'), ylabel('y')


%% 
%parameters
t_end = 30;
options = odeset('RelTol',1e-8,'Maxstep',0.1);
amplitudes = 0:.2:1;
for a=1:length(amplitudes)
    [t,y] = ode15s(@(t,y) FHNExcitability(t,y,amplitudes(a)),[0 t_end],y0,options);
    %plot numerical solution
    figure(1)
    subplot(2,1,1);
    plot(t,y(:,1),'-');
    hold on;
    ylabel('x');
    set(gca,'Fontsize',14,'LineWidth',1)
    subplot(2,1,2);
    plot(t,y(:,2),'-');
    hold on;
    
    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('Time (t.u.)'), ylabel('y')
    title('Excitability')
    
    %plot numerical solution
    figure(2)
    plot(y(:,1),y(:,2),'-','LineWidth',2)
    hold on;
    
    x_range = -2.5:.1:2;
    nullcline_1 = (1/0.8)*(x_range+0.7);
    nullcline_2 = (-1/3)*x_range.^3+x_range;
    plot(x_range, nullcline_1, 'g--')
    plot(x_range, nullcline_2, 'k--')
    title('Excitability')
    hold on
    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('x'), ylabel('y')

end
figure(1);
subplot(2,1,1);
legend(strcat('amp=',string(num2cell(amplitudes))))

%% 
t_end = 50;
options = odeset('RelTol',1e-8,'Maxstep',0.1);
interval = 10;

[t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
%plot numerical solution
figure(1)
subplot(2,1,1);
plot(t,y(:,1),'-');
hold on;
plot(5, -1.9, 'rv','LineWidth',1);
plot(5+interval, -1.9, 'rv','LineWidth',1);
title('Successive Stimulations')
legend('', 'Stimulation')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
subplot(2,1,2);
plot(t,y(:,2),'-');
hold on;
plot( 5, -.95, 'rv','LineWidth',1);
plot( 5+interval, -.95, 'rv','LineWidth',1);
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')

%% 

t_end = 50;
options = odeset('RelTol',1e-8,'Maxstep',0.1);
interval = 9;

[t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
%plot numerical solution
figure(1)
subplot(2,2,1);
plot(t,y(:,1),'-');
hold on;
plot(5, -1.9, 'rv','LineWidth',1);
plot(5+interval, -1.9, 'rv','LineWidth',1);
title('Interval=9')
legend('', 'Stimulation')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
subplot(2,2,3);
plot(t,y(:,2),'-');
hold on;
plot( 5, -.95, 'rv','LineWidth',1);
plot( 5+interval, -.95, 'rv','LineWidth',1);
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')

interval=8;
[t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
%plot numerical solution
figure(1)
subplot(2,2,2);
plot(t,y(:,1),'-');
hold on;
plot(5, -1.9, 'rv','LineWidth',1);
plot(5+interval, -1.9, 'rv','LineWidth',1);
title('Interval=8')
legend('', 'Stimulation')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
subplot(2,2,4);
plot(t,y(:,2),'-');
hold on;
plot( 5, -.95, 'rv','LineWidth',1);
plot( 5+interval, -.95, 'rv','LineWidth',1);
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')

%% 
t_end = 30;
interval = 9;
[t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
figure(1)
subplot(2,2,1);
plot(t,y(:,1),'-');
hold on;
plot(5, -1.9, 'rv','LineWidth',1);
plot(5+interval, -1.9, 'rv','LineWidth',1);
title('Interval=9')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
interval = 10;
[t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
figure(1)
subplot(2,2,2);
plot(t,y(:,1),'-');
hold on;
plot(5, -1.9, 'rv','LineWidth',1);
plot(5+interval, -1.9, 'rv','LineWidth',1);
title('Interval=10')
legend('', 'Stimulation')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)

interval = 11;
[t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
figure(1)
subplot(2,2,3);
plot(t,y(:,1),'-');
hold on;
plot(5, -1.9, 'rv','LineWidth',1);
plot(5+interval, -1.9, 'rv','LineWidth',1);
title('Interval=11')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)

interval = 12;
[t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
figure(1)
subplot(2,2,4);
plot(t,y(:,1),'-');
hold on;
plot(5, -1.9, 'rv','LineWidth',1);
plot(5+interval, -1.9, 'rv','LineWidth',1);
title('Interval=12')
hold on;
ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)

%% 
t_end = 100;
interval_values = 9.01:14;
latencies = zeros(length(interval_values),1);
for v=1:length(interval_values)
    interval = interval_values(v);
    [t,y] = ode15s(@(t,y) FHNRefrac(t,y, interval),[0 t_end],y0,options);
    [pks, locs] = findpeaks(y(:,1));
    peak_time = t(locs(length(pks)));
    latency = peak_time - (5+interval);
    latencies(v) = latency;
end
figure();
plot(interval_values, latencies);
xlabel('Interval between stims (t.u)')
ylabel('Time from stim to peak (t.u.)')
title('Latency')
set(gca,'Fontsize',14,'LineWidth',1)

%% 


% initial conditions
y0 = [-1.2; -0.62];

%parameters
t_end = 50;

options = odeset('RelTol',1e-8,'Maxstep',0.1);
[t,y] = ode15s(@(t,y) FHNExcitability(t,y,-2.5),[0 t_end],y0,options);
%plot numerical solution
figure(1)
subplot(2,1,1);
plot(t,y(:,1),'b-');
ylim([-2.5 2])
hold on;
plot(5, -2.4, 'r*','LineWidth',1);
legend('', 'Hyperpolarizing Stim (Amp = -2.5)')
title('Anodal Break Response')

ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
subplot(2,1,2);
plot(t,y(:,2),'b-');
ylim([-1.2 1.5])
hold on;
plot(5, -1.15, 'r*','LineWidth',1);

set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')


%plot numerical solution
figure(2)
plot(y(:,1),y(:,2),'b-','LineWidth',2)
hold on;

x_range = -2.5:.1:2;
nullcline_1 = (1/0.8)*(x_range+0.7);
nullcline_2 = (-1/3)*x_range.^3+x_range;
plot(x_range, nullcline_1, 'g--')
plot(x_range, nullcline_2, 'k--')

legend('', 'Ydot nullcline', 'Xdot nullcline')
title('Hyperpolarizing Current')
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('x'), ylabel('y')
%% % initial conditions
y0 = [-1.2; -0.62];

%parameters
t_end = 100;

a = 0;

options = odeset('RelTol',1e-8,'Maxstep',0.1);
[t,y] = ode15s(@(t,y) FHNA(t,y, a),[0 t_end],y0,options);
%plot numerical solution
figure(1)
subplot(2,1,1);
plot(t,y(:,1),'b-');
hold on;
title('a=0')

ylabel('x');
set(gca,'Fontsize',14,'LineWidth',1)
subplot(2,1,2);
plot(t,y(:,2),'b-');
hold on;

set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)'), ylabel('y')


%plot numerical solution
figure(2)
plot(y(:,1),y(:,2),'b-','LineWidth',2)
hold on;

x_range = -2.5:.1:2;
nullcline_1 = (1/0.8)*(x_range+a);
nullcline_2 = (-1/3)*x_range.^3+x_range;
% plot(x_range, nullcline_1, 'g--')
% plot(x_range, nullcline_2, 'k--')

% legend('', 'Ydot nullcline', 'Xdot nullcline')
title('a=0')
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('x'), ylabel('y')

%% 
y0 = [-1.2; -0.62];

%parameters
t_end = 100;

a_values = [.5, .01, -.01, -.5];

options = odeset('RelTol',1e-8,'Maxstep',0.1);
for a_idx=1:length(a_values)
    [t,y] = ode15s(@(t,y) FHNA(t,y, a_values(a_idx)),[0 t_end],y0,options);
    %plot numerical solution
    figure(1)
    subplot(2,1,1);
    plot(t,y(:,1),'-');
    hold on;
    
    ylabel('x');
    set(gca,'Fontsize',14,'LineWidth',1)
    subplot(2,1,2);
    plot(t,y(:,2),'-');
    hold on;
    
    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('Time (t.u.)'), ylabel('y')
    
    
    %plot numerical solution
    figure(2)
    plot(y(:,1),y(:,2),'-','LineWidth',1)
    hold on;

    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('x'), ylabel('y')
end

figure(1);
subplot(2,1,1);
legend(strcat('a=',string(num2cell(a_values))));
figure(2);
legend(strcat('a=',string(num2cell(a_values))));

%% 
y0 = [-1.2; -0.62];

%parameters
t_end = 150;

a_values = [.44, .43, .42, .41];

options = odeset('RelTol',1e-8,'Maxstep',0.1);
for a_idx=1:length(a_values)
    [t,y] = ode15s(@(t,y) FHNA(t,y, a_values(a_idx)),[0 t_end],y0,options);
    %plot numerical solution
    figure(1)
    subplot(2,1,1);
    plot(t,y(:,1),'-');
    hold on;
    
    ylabel('x');
    set(gca,'Fontsize',14,'LineWidth',1)
    subplot(2,1,2);
    plot(t,y(:,2),'-');
    hold on;
    
    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('Time (t.u.)'), ylabel('y')
    
    
    %plot numerical solution
    figure(2)
    plot(y(:,1),y(:,2),'-','LineWidth',1)
    hold on;

    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('x'), ylabel('y')
end

figure(1);
subplot(2,1,1);
legend(strcat('a=',string(num2cell(a_values))));
figure(2);
legend(strcat('a=',string(num2cell(a_values))));

%% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = FHNeqns(t,y)

a = 0.7; 
b = 0.8;
c = 3.0;

amp = .0;
dur = 0.5;
tstart = 5.;
S=0;
if (t>tstart && t<tstart+dur) 
    S=amp;
end

dydt = [c*(y(1)-1/3*y(1)^3-y(2)+S); (y(1)+a-b*y(2))/c];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = FHNExcitability(t,y, amp)

a = 0.7; 
b = 0.8;
c = 3.0;


dur = 0.75;
tstart = 5.;
S=0;
if (t>tstart && t<tstart+dur) 
    S=amp;
end

dydt = [c*(y(1)-1/3*y(1)^3-y(2)+S); (y(1)+a-b*y(2))/c];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = FHNRefrac(t,y, interval)

a = 0.7; 
b = 0.8;
c = 3.0;
amp=0.7;

dur = 0.5;
tstart = 5.;
tstart2 = tstart + interval;
S=0;
if (t>tstart && t<tstart+dur) 
    S=amp;
end
if (t>tstart2 && t<tstart2+dur)
    S=amp;
end

dydt = [c*(y(1)-1/3*y(1)^3-y(2)+S); (y(1)+a-b*y(2))/c];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = FHNA(t,y, a)

b = 0.8;
c = 3.0;

amp = .0;
dur = 0.5;
tstart = 5.;
S=0;
if (t>tstart && t<tstart+dur) 
    S=amp;
end

dydt = [c*(y(1)-1/3*y(1)^3-y(2)+S); (y(1)+a-b*y(2))/c];

end
