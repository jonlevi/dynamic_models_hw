options = odeset('RelTol',1e-8);

t_end = 200;

params = {};
params.vs = 0.76;
params.kI = 1;
params.n = 4;
params.vm = 0.65;
params.km1 = 0.5;
params.ks = 0.38;
params.v1 = 3.2;
params.k1 = 1;
params.k2 = 1;
params.v2 = 1.58;
params.v3 = 5;
params.v4 = 2.5;
params.k3 = 2;
params.k4 = 2;
params.kd = 0.2;
params.vd = 0.95;

params.kp1 = 1.9;
params.kp2 = 1.3;

y0 = [1 1 1 1 1];

[t,y] = ode45(@(t,y) model(t,y, params),[0 t_end],y0,options);

figure();

for i=1:5
    plot(t, y(:,i));
    hold on;
end
xlabel('time (hr)')
ylabel('Concentration (\muM)')
legend('per mRNA', 'p0', 'p1', 'p2', 'pN');

figure();
plot(t, y(:,1));
hold on;
plot(t, y(:,2) + y(:,3) + y(:,4) + y(:,5));
plot(t, y(:,5));
xlabel('time (hr)')
ylabel('Concentration (\muM)')
legend('per mRNA', 'Total PER', 'Nuclear PER');

figure();
plot(t, y(:,1));
[pks,locs] = findpeaks(y(:,1));
xlabel('time (hr)')
ylabel('mRNA Concentration (\muM)')
hold on;
plot(t(locs), pks, 'o');
t(locs)
plot([t(locs(6)) t(locs(7))], [pks(6) pks(7)], 'k--')
xlim([t(locs(6))-10 , t(locs(7))+10]);

figure();
plot(t, y(:,2));
hold on;
plot(t, y(:,4));
xlabel('time (hr)')
ylabel('Concentration (\muM)')
legend('p0', 'p2');
%% 
vd_values = .45:.1:2.6;
period_values = zeros(length(vd_values),1);
figure();
for i=1:length(period_values)
    params.vd = vd_values(i);
    [t,y] = ode45(@(t,y) model(t,y, params),[0 t_end],y0,options);
    [pks,locs] = findpeaks(y(:,1));
    period_values(i) = mean(diff(t(locs)));

end

params.vd = 0.95;

figure(1);
plot(vd_values, period_values);
xlabel('v_d')
ylabel('Period of per mRNA oscillation')

%% 
vd_values = 2.53:.01:2.56;
amp_values = zeros(length(vd_values),1);
for i=1:length(amp_values)
    params.vd = vd_values(i);
    [t,y] = ode45(@(t,y) model(t,y, params),[0 t_end*5],y0,options);
    plot(t, y(:,1));
    hold on;
    [pks,locs] = findpeaks(y(:,1));
    amp_values(i) = mean(pks);

end
xlabel('time (hr)')
ylabel('per mRNA concentration')
legend(strcat('v_d=',string(num2cell(vd_values))))
params.vd = 0.95;

figure();
% plot(sqrt(2.6-vd_values), amp_values);
plot((vd_values), amp_values);
xlabel('sqrt(2.6-v_d)')
ylabel('Amplitude of per mRNA oscillation')

%% 

vd_values = .47:.001:.48;
amp_values = zeros(length(vd_values),1);
for i=1:length(amp_values)
    params.vd = vd_values(i);
    [t,y] = ode45(@(t,y) model(t,y, params),[0 t_end*100],y0,options);
    plot(t, y(:,1));
    hold on;
    transient = y(t>500,1);
    [pks,locs] = findpeaks(transient);
    [pks2,locs2] = findpeaks(-1*transient);
    amp_values(i) = pks(end)+pks2(end);

end
xlabel('time (hr)')
ylabel('per mRNA concentration')
legend(strcat('v_d=',string(num2cell(vd_values))))
params.vd = 0.95;

figure();
% plot((vd_values), log10(amp_values));
plot(sqrt(vd_values), amp_values);
xlabel('sqrt(v_d)')
ylabel('Amplitude of per mRNA oscillation')


%% 
n_values = [3 4 6 8];
vd_values = .45:.1:2.6;
period_values = zeros(length(n_values), length(vd_values));
for i=1:length(n_values)
    
    params.n = n_values(i);
    % [t,y] = ode45(@(t,y) model(t,y, params),[0 t_end*2],y0,options);
    % plot(t, y(:,1));
    % hold on;
    
    % period_values = zeros(length(vd_values),1);
    % figure();
    for j=1:length(vd_values)
        params.vd = vd_values(j);
        [t,y] = ode45(@(t,y) model(t,y, params),[0 2*t_end],y0,options);
        % figure(1+i);
        % plot(t, y(:,1));
        hold on;
        [pks,locs] = findpeaks(y(:,1));
        period_values(i,j) = mean(diff(t(locs)));
    
    end
    % legend(strcat('vd=',string(num2cell(vd_values))))

    params.vd = 0.95;


end
% legend(strcat('n=',string(num2cell(n_values))))
% xlabel('time (hr)')
% ylabel('per mRNA concentration')
params.n = 4;
% figure();
% plot(n_values, period_values);
% xlabel('n')
% ylabel('Period of per mRNA oscillation')
figure();
for x=1:length(n_values)
    plot(vd_values, period_values(x,:));
    hold on;
end
legend(strcat('n=',string(num2cell(n_values))))
xlabel('v_d')
ylabel('Period of per mRNA oscillation')

%% 
y0 = [1 1 1];

[t,y] = ode45(@(t,y) one_per_model(t,y, params),[0 t_end],y0,options);

figure();

for i=1:3
    plot(t, y(:,i));
    hold on;
end
xlabel('time (hr)')
ylabel('Concentration (\muM)')
legend('per mRNA', 'pC', 'pN');



%% 
y0 = [1 1 1];
vals = [0.01];
for blah = 1:length(vals)
    params.kd = vals(blah);
    [t,y] = ode45(@(t,y) one_per_model(t,y, params),[0 t_end],y0,options);
    params.kd = 0.2;
    figure();
    
    for i=1:3
        plot(t, y(:,i));
        hold on;
    end
    xlabel('time (hr)')
    ylabel('Concentration (\muM)')
    legend('per mRNA', 'pC', 'pN');

end


figure();
pc = 0:.1:1;
params.kd = 0.01;
plot(pc, (params.vd*pc)./(params.kd+pc));
hold on;
params.kd = 0.2;
plot(pc, (params.vd*pc)./(params.kd+pc));
ylabel('Degradation Term: (V_d p_C) / (K_d + p_C)')
xlabel('p_C')
legend('K_d = 0.01', 'K_d = 0.2')

%% 
lag = 4;
figure();
y0 = [1 1 1];
opts = ddeset('RelTol',1e-8,'AbsTol',1e-8);
sol =  dde23(@(t,y,Z) delay_model(t,y,Z,params),lag,@(t) y0,[0, t_end],opts);
for i=1:3
    plot(sol.x, sol.y(i,:));
    hold on;
end
xlabel('time (hr)')
ylabel('Concentration (\muM)')
legend('per mRNA', 'pC', 'pN');

%% 
lag = 4;
y0 = [1 1 1];
opts = ddeset('RelTol',1e-8,'AbsTol',1e-8);


vd_values = .47:.1:2.5;
period_values = zeros(length(vd_values),1);

for i=1:length(period_values)
    % figure();
    params.vd = vd_values(i);
    sol =  dde23(@(t,y,Z) delay_model(t,y,Z,params),lag,@(t) y0,[0, t_end],opts);
    [pks,locs] = findpeaks(sol.y(1,:));
    period_values(i) = mean(diff(sol.x(locs)));
    for j=1:3
        % plot(sol.x, sol.y(j,:));
        hold on;
    end
% title(strcat('vd=', num2str(params.vd)))
% xlabel('time (hr)')
% ylabel('Concentration (\muM)')
% legend('per mRNA', 'pC', 'pN');

end

params.vd = 0.95;

figure(1);
hold on;
plot(vd_values, period_values);
xlabel('v_d')
ylabel('Period of per mRNA oscillation')
%% 


y0 = [1 1 1];
opts = ddeset('RelTol',1e-8,'AbsTol',1e-8);
figure();
lags = 1.4:.05:1.6;
for li=1:length(lags)
    sol =  dde23(@(t,y,Z) delay_model(t,y,Z,params),lags(li),@(t) y0,[0, 5*t_end],opts);
    plot(sol.x, sol.y(1,:));
    hold on;
end
legend(strcat('\tau=',string(num2cell(lags))))

xlabel('time (hr)')
ylabel('Concentration of per mRNA (\muM)')


%% 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = model(t,y, params)

m = y(1);
p0 = y(2);
p1 = y(3);
p2 = y(4);
pn = y(5);

v1_term = (params.v1*p0)/(params.k1+p0);
v2_term = (params.v2*p1)/(params.k2+p1);
v3_term = (params.v3*p1)/(params.k3+p1);
v4_term = (params.v4*p2)/(params.k4+p2);
vd_term = (params.vd*p2)/(params.kd+p2);
vm_term = (params.vm*m)/(params.km1 + m);
vs_term = params.vs/((1+ (pn/params.kI).^params.n));

dydt = [vs_term - vm_term;
       params.ks*m - v1_term + v2_term;
        v1_term - v2_term - v3_term + v4_term;
        v3_term - v4_term - vd_term - params.kp1*p2 + params.kp2*pn;
        params.kp1*p2 - params.kp2*pn];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = one_per_model(t,y, params)

m = y(1);
pc = y(2);
pn = y(3);

vd_term = (params.vd*pc)/(params.kd+pc);
vm_term = (params.vm*m)/(params.km1 + m);
vs_term = params.vs/((1+ (pn/params.kI).^params.n));

dydt = [vs_term - vm_term;
       params.ks*m - vd_term - params.kp1*pc + params.kp2*pn;
        params.kp1*pc - params.kp2*pn];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dydt = delay_model(t,y,Z, params)
m = y(1);
pc = y(2);
pn = y(3);

ylag = Z(:,1);
pc_delay = ylag(2);

vd_term = (params.vd*pc)/(params.kd+pc);
vm_term = (params.vm*m)/(params.km1 + m);
vs_term = params.vs/((1+ (pn/params.kI).^params.n));


dydt = [vs_term - vm_term;
       params.ks*m - vd_term - params.kp1*pc_delay + params.kp2*pn;
        params.kp1*pc_delay - params.kp2*pn];

end