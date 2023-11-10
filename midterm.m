options = odeset('RelTol',1e-8, 'NonNegative',1);

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;


initial_values = [1.5*10^4 1000 1000 1000];
initial_values2 = [1.5*10^2 1000 1000 1000];
initial_values3 = [1.5*10^3 1000 1000 1000];
initial_values4 = [1.5*10^4 1000 1000 1000];
[t,y] = ode45(@(t,y) model(t,y,params),[0 200000],initial_values,options);
[t2,y2] = ode45(@(t,y) model(t,y,params),[0 200000],initial_values2,options);
[t3,y3] = ode45(@(t,y) model(t,y,params),[0 200000],initial_values3,options);
[t4,y4] = ode45(@(t,y) model(t,y,params),[0 200000],initial_values4,options);

t = t/10^4;
t2 = t2/10^4;
t3 = t3/10^4;
t4 =  t4/10^4;




figure();
subplot(4,1,1);
plot(t, y(:,1), 'r');
hold on;
% plot(t2, y2(:,1), 'b--');
hold on;
% plot(t3, y3(:,1), 'g--');
% hold on;
% plot(t4, y4(:,1), 'c--');
% ylim([-50000 5*10^5]);
ylabel('Mosquito Population (M)');
% legend('Low M_0', 'Medium M_0', 'High M_0', 'Very High M_0');
subplot(4,1,2);
plot(t, y(:,2), 'r');
hold on;
% plot(t2, y2(:,2), 'b--');
% plot(t3, y3(:,2), 'g--');
% plot(t4, y4(:,2), 'c--');
% ylim([-20 400]);

ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'r');
% hold on;
% plot(t2, y2(:,3), 'b--');
% plot(t3, y3(:,3), 'g--');
% plot(t4, y4(:,3), 'c--');
% ylim([-5000 5*10^4])

ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'r');
hold on;
% plot(t2, y2(:,4), 'b--');
% plot(t3, y3(:,4), 'g--');
% plot(t4, y4(:,4), 'c--');
% ylim([200 700]);

ylabel('Rat Population (R)');
% xlim([49.9 50])
xlabel('Time (a.u.)');
% sgtitle('Stable Attractors (parameter set #1)')


%% 

params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;


parameter_scaling = 0.1:0.1:2;
initial_values = [1.5*10^4 1000 1000 1000];
param_values = zeros(length(parameter_scaling),1);
final_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.rm = 0.35*parameter_scaling(i);
    param_values(i) = 0.35*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
figure(1);
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');

params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;


param_values = zeros(length(parameter_scaling),1);
final_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.km = 2*10^5*parameter_scaling(i);
    param_values(i) = 0.35*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value');
ylabel('Steady State Value of M');

params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;

param_values = zeros(length(parameter_scaling),1);
final_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.am = 0.35*parameter_scaling(i);
    param_values(i) = 0.35*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end

plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value');
ylabel('Steady State Value of M');

params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;

param_values = zeros(length(parameter_scaling),1);
final_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.bm = 0.5*parameter_scaling(i);
    param_values(i) = 0.35*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end

plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value');
ylabel('Steady State Value of M');

params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;

param_values = zeros(length(parameter_scaling),1);
final_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.dm = 0.05*parameter_scaling(i);
    param_values(i) = 0.35*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value');
ylabel('Steady State Value of M');

%% 
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;


parameter_scaling = 0.1:0.1:2;
initial_values = [1.5*10^4 1000 1000 1000];
final_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.qs = 0.75*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
figure(2);
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');

final_values = zeros(length(parameter_scaling),1);
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;

for i=1:length(parameter_scaling)
    params.ds = 0.001*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');
%

final_values = zeros(length(parameter_scaling),1);
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;
for i=1:length(parameter_scaling)
    params.alpha = .00002*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');

%
final_values = zeros(length(parameter_scaling),1);
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;
for i=1:length(parameter_scaling)
    params.beta = .000002*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');

%
final_values = zeros(length(parameter_scaling),1);
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;
for i=1:length(parameter_scaling)
    params.dl =  0.008*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');
%% 
figure(2);
parameter_scaling = 0.5:0.2:1.5;
final_values = zeros(length(parameter_scaling),1);
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;
for i=1:length(parameter_scaling)
    params.rr =  0.25*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');

final_values = zeros(length(parameter_scaling),1);
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;
for i=1:length(parameter_scaling)
    params.kr =  500*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');

final_values = zeros(length(parameter_scaling),1);
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;
for i=1:length(parameter_scaling)
    params.dr =  .05*parameter_scaling(i);
    [~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
    final_values(i) = y(length(y));

end
plot(parameter_scaling, final_values);
hold on;
xlabel('Parameter Value')
ylabel('Steady State Value of M');
%% 
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;

param_base = 0.05;
param_plus = param_base + .1*param_base;
param_minus = param_base - .1*param_base;

[~,y] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
mbase = y(length(y));
params.dr = param_plus;
[~,yplus] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
mplus = yplus(length(yplus));
params.dr = param_minus;
[~,yminus] = ode45(@(t,y) model(t,y,params),[0 500000],initial_values,options);
mminus = yminus(length(yminus));

dp = param_plus - param_minus;
ds = mplus - mminus;
p = param_base;
s = mbase;
disp('sensitivity:')
disp((p*ds)/(dp*s));

%% 


%% 
initial_values = [1.5*10^5 1000 1000 1000];
[t1,y1] = ode45(@(t,y) modelDDT(t,y,params, -14000),[0 500000],initial_values,options);
plot(t1, y1(:,1), 'k--');
hold on;
[t,y] = ode45(@(t,y) modelDDT(t,y,params, -12500),[0 500000],initial_values,options);
plot(max(t1):max(t), zeros(length((max(t1):max(t))),1), 'k--')
xlim([0 1.5*10^5])
plot(t, y(:,1), 'r');
ylim([-50000 5*10^5]);
ylabel('Mosquito Population (M)');
% subplot(4,1,2);
% plot(t, y(:,2), '--');
% hold on;
% 
% ylabel('Small Bird Population (S)');
% 
% subplot(4,1,3);
% plot(t, y(:,3),'--');
% hold on;
% 
% ylabel('Large Bird Population (L)');
% 
% subplot(4,1,4);
% plot(t, y(:,4), '--');
% hold on;
% 
% ylabel('Rat Population (R)');
% xlabel('Time (a.u.)');

%% 
options = odeset('RelTol',1e-8, 'NonNegative',1);

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;


initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
[p1,f1] = pspectrum(prey(:,1),1);
pret = pret/10^4;
cutoff = ceil(length(pret)/500);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'r');
hold on;
ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'r');
hold on;
ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);

ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'r');
hold on;
% ylim([900 3000]);
ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'r');
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');



params.dr = 0.18;


initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
pret = pret/10^4;
cutoff = ceil(length(pret)/500);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'b');
hold on;
ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'b');
hold on;
ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);

ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'b');
% ylim([900 3000]);
ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'b');
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');

params.dr = 0.185;


initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
[p2,f2] = pspectrum(prey(:,1),1);
pret = pret/10^4;
cutoff = ceil(length(pret)/500);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'k');
legend('d_R = 0.16', 'd_R = 0.18', 'd_R = 0.185')
hold on;
ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'k');
hold on;
ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);
xlim([min(t) max(t)])
ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'k');
% ylim([900 3000]);
ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
xlim([min(t) max(t)])
ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'k');
xlim([min(t) max(t)])
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');

% figure();
% plot(1/length(t)*(0:length(t)-1),abs(fft(y(:,1))),"LineWidth",1);
% hold on;
% plot(1/length(t)*(0:length(t)-1),abs(fft(y(:,2))),"LineWidth",1)
% plot(1/length(t)*(0:length(t)-1),abs(fft(y(:,3))),"LineWidth",1)
% 
% plot(1/length(t)*(0:length(t)-1),abs(fft(y(:,4))),"LineWidth",1)
%% 
options = odeset('RelTol',1e-8, 'NonNegative',1);

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.05;

params.dm = 0.05;


initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
pret = pret/10^4;
cutoff = ceil(length(pret)/100);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'r', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'r', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);

ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'r', 'MarkerSize', .5);
hold on;
% ylim([900 3000]);
% ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
% % ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'r', 'MarkerSize', .5);
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');


params.dm = 0.08;
params.ds = 0.0015;
params.dl = 0.01;


initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
pret = pret/10^4;
cutoff = ceil(length(pret)/100);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'b--', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'b--', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);

ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'b--', 'MarkerSize', .5);
% ylim([900 3000]);
% ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'b--', 'MarkerSize', .5);
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');

params.dm = 0.1;
params.ds = 0.002;
params.dl = 0.015;



initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
pret = pret/10^4;
cutoff = ceil(length(pret)/100);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'k--');
% legend('d_m = 0.05', 'd_m = 0.1', 'd_m = 0.3')
hold on;
% ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'k--');
hold on;
% ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);
% xlim([min(t) max(t)])
ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'k--');
% ylim([900 3000]);
% ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
% xlim([min(t) max(t)])
ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'k--');
% xlim([min(t) max(t)])
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');

%% 
%% 
options = odeset('RelTol',1e-8, 'NonNegative',1);

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;

params.dm = 0.05;


initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
pret = pret/10^4;
cutoff = ceil(length(pret)/100);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'r', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'r', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);

ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'r', 'MarkerSize', .5);
hold on;
% ylim([900 3000]);
% ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
% % ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'r', 'MarkerSize', .5);
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');


params.dm = 0.08;
params.ds = 0.0015;
params.dl = 0.01;


initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
pret = pret/10^4;
cutoff = ceil(length(pret)/100);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'b--', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'b--', 'MarkerSize', .5);
hold on;
% ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);

ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'b--', 'MarkerSize', .5);
% ylim([900 3000]);
% ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'b--', 'MarkerSize', .5);
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');

params.dm = 0.1;
params.ds = 0.002;
params.dl = 0.015;



initial_values = [1.5*10^4 1000 1000 1000];
[pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
pret = pret/10^4;
cutoff = ceil(length(pret)/100);
t = pret(length(pret)-cutoff:length(pret));
y = prey(length(pret)-cutoff:length(pret), :);
figure(4);
subplot(4,1,1);
plot(t, y(:,1), 'k--');
% legend('d_m = 0.05', 'd_m = 0.1', 'd_m = 0.3')
hold on;
% ylim([min(y(:,1))*.999 max(y(:,1))*1.001]);
ylabel('Mosquito Population (M)');
subplot(4,1,2);
plot(t, y(:,2), 'k--');
hold on;
% ylim([min(y(:,2))*.9 max(y(:,2))*1.1]);
% xlim([min(t) max(t)])
ylabel('Small Bird Population (S)');

subplot(4,1,3);
plot(t, y(:,3),'k--');
% ylim([900 3000]);
% ylim([min(y(:,3))*.9 max(y(:,3))*1.1]);
% xlim([min(t) max(t)])
ylabel('Large Bird Population (L)');

subplot(4,1,4);
plot(t, y(:,4), 'k--');
% xlim([min(t) max(t)])
% ylim([10 50]);
% ylim([min(y(:,4))*.9 max(y(:,4))*1.1]);
hold on;

ylabel('Rat Population (R)');
xlabel('Time (a.u.)');


%% 

options = odeset('RelTol',1e-8, 'NonNegative',1);
initial_values = [1.5*10^4 1000 1000 1000];

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;

figure()
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.rm =  .35*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end

legend(strcat('SF=',string(num2cell(parameter_scaling))))
title('changing rm')

figure(2);
hold on;
plot(parameter_scaling, amp_values);
hold on;
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
hold on;
plot(parameter_scaling,period_values);
hold on;
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')


params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.km =  2*10^5*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing km')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.am =  0.35*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing am')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')


params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.bm =  0.5*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing bm')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.dm = 0.05*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing dm')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.qs =  0.75*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing qs')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
hold on;
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
hold on;
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.ds = 0.001*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing ds')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')



params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.alpha = .00002*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing alpha')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.beta = .000002*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing beta')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')

params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    params.dl = 0.008*parameter_scaling(i);
    [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
    cutoff = ceil(length(pret)/1000);
    t = pret(length(pret)-cutoff:length(pret));
    y = prey(length(pret)-cutoff:length(pret), :);
    plot(t, y(:,4));
    hold on;
    ylabel('Rat Population (R)');
    xlabel('Time (a.u.)');
    [pks,locs] = findpeaks(y(:,4));
    [minpks,minlocs] = findpeaks(-1*y(:,4));
    amp_values(i) = mean(pks) - mean(-1*minpks);
    period_values(i) = mean(diff(t(locs)));
end
title('changing dl')
legend(strcat('SF=',string(num2cell(parameter_scaling))))

figure(2);
plot(parameter_scaling, amp_values);
xlabel('Parameter Scaling')
ylabel('Amplitude of R Oscillation')
figure(3);
plot(parameter_scaling,period_values);
xlabel('Parameter Scaling')
ylabel('Period of R Oscillation')



%% 
params = {};
params.rm = 0.35;
params.km = 2*10^5;
params.am = 0.35;
params.bm = 0.5;
params.dm = 0.05;
params.qs = 0.75;
params.ds = 0.001;
params.alpha = .00002;
params.beta = .000002;
params.dl = 0.008;
params.rr = 0.25;
params.kr = 500;
params.dr = 0.18;



figure();
amp_values = zeros(length(parameter_scaling),1);
period_values = zeros(length(parameter_scaling),1);

for i=1:length(parameter_scaling)
    for j=1:length(parameter_scaling)
        params.qs = 0.75*parameter_scaling(i);
        params.beta = .000002*parameter_scaling(j);
        [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
        cutoff = ceil(length(pret)/1000);
        t = pret(length(pret)-cutoff:length(pret));
        y = prey(length(pret)-cutoff:length(pret), :);
        plot(t, y(:,4));
        hold on;
        ylabel('Rat Population (R)');
        xlabel('Time (a.u.)');
        [pks,locs] = findpeaks(y(:,4));
        [minpks,minlocs] = findpeaks(-1*y(:,4));
        amp_values(i,j) = mean(pks) - mean(-1*minpks);
        period_values(i,j) = mean(diff(t(locs)));
    end
end

figure(2);
surf(parameter_scaling,parameter_scaling,amp_values);
colorbar;
xlabel('q_S Parameter Scaling')
ylabel('\beta Parameter Scaling')
zlabel('Amplitude of R Oscillation')

figure(3);
colorbar;
surf(parameter_scaling,parameter_scaling,period_values);
xlabel('q_S Parameter Scaling')
ylabel('\beta Parameter Scaling')
zlabel('Period of R Oscillation')



% figure();
% amp_values = zeros(length(parameter_scaling),1);
% period_values = zeros(length(parameter_scaling),1);
% 
% for i=1:length(parameter_scaling)
%     params.beta = .000002*parameter_scaling(i);
%     [pret,prey] = ode45(@(t,y) model(t,y,params),linspace(0, 500000, 500000),initial_values,options);
%     cutoff = ceil(length(pret)/1000);
%     t = pret(length(pret)-cutoff:length(pret));
%     y = prey(length(pret)-cutoff:length(pret), :);
%     plot(t, y(:,4));
%     hold on;
%     ylabel('Rat Population (R)');
%     xlabel('Time (a.u.)');
%     [pks,locs] = findpeaks(y(:,4));
%     [minpks,minlocs] = findpeaks(-1*y(:,4));
%     amp_values(i) = mean(pks) - mean(-1*minpks);
%     period_values(i) = mean(diff(t(locs)));
% end
% title('changing beta')
% legend(strcat('SF=',string(num2cell(parameter_scaling))))
% 
% figure(2);
% plot(parameter_scaling, amp_values);
% xlabel('Parameter Scaling')
% ylabel('Amplitude of R Oscillation')
% figure(3);
% plot(parameter_scaling,period_values);
% xlabel('Parameter Scaling')
% ylabel('Period of R Oscillation')
%% 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = model(t,y, params)

M = y(1);
S = y(2);
L = y(3);
R = y(4);

m_prey_term = (params.am*M)/(1+params.bm*M);

dydt = [params.rm * M*(1-(M/params.km)) - S*m_prey_term - params.dm*M;
        S*(params.qs * m_prey_term - (params.alpha*L) - params.ds);
        L*(params.alpha*S + params.beta*R - params.dl)
        params.rr * R * (1-(R/params.kr)) - params.beta*R*L - params.dr*R
        ];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = modelDDT(t,y, params, ddt)


M = max(0, y(1));
dec=0;
if (t > 10*10^4) && ((t < 10.01*10^4))
    dec = ddt;
end
S = y(2);
L = y(3);
R = y(4);

m_prey_term = (params.am*M)/(1+params.bm*M);

dydt = [dec+(params.rm * M*(1-(M/params.km)) - S*m_prey_term - params.dm*M);
        S*(params.qs * m_prey_term - (params.alpha*L) - params.ds);
        L*(params.alpha*S + params.beta*R - params.dl)
        params.rr * R * (1-(R/params.kr)) - params.beta*R.*L - params.dr*R
        ];

end
