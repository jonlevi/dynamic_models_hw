


k1_def = 0.15;
%k1 = 0.15;
k2 = 0.2;
k3 = 0.25;


% initial conditions
y0 = [0; 1; 1];

%parameters
t_end = 5000;
t_transient = 4500;

sensfac = 0.1:0.1:2;

maxG = zeros(length(sensfac),1);
minG = zeros(length(sensfac),1);

for i =1:length(sensfac)
    
    k1 = sensfac(i)*k1_def;

    options = odeset('RelTol',1e-8);
    [t,y] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);

    %determine max and min of G, after transient
    mini = min(find(t > t_transient)); 
    l = length(t);
    maxG(i) = max(y(mini:l,3));
    minG(i) = min(y(mini:l,3));
        
end

%plot numerical solution
figure()
subplot(2,1,1);
plot(sensfac,maxG,'b*-');
hold on;
plot(sensfac,minG,'r*-');
legend('MaxG', 'MinG')
set(gca,'Fontsize',14,'LineWidth',1)
ylabel('G max and min')


subplot(2,1,2);
plot(sensfac,maxG-minG,'b*-')
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('k1 parameter scaling')
ylabel('G amplitude (maxG-minG)')
sgtitle('K1 Sensitivity')

%% 
k1 = 0.15;
k2_def = 0.2;
k3 = 0.25;
% initial conditions
y0 = [0; 1; 1];

%parameters
t_end = 5000;
t_transient = 4500;

sensfac = 0.1:0.1:2;

maxG = zeros(length(sensfac),1);
minG = zeros(length(sensfac),1);

for i =1:length(sensfac)
    
    k2 = sensfac(i)*k2_def;

    options = odeset('RelTol',1e-8);
    [t,y] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);

    %determine max and min of G, after transient
    mini = min(find(t > t_transient)); 
    l = length(t);
    maxG(i) = max(y(mini:l,3));
    minG(i) = min(y(mini:l,3));
        
end

%plot numerical solution
figure()
subplot(2,1,1);
plot(sensfac,maxG,'b*-');
hold on;
plot(sensfac,minG,'r*-');
legend('MaxG', 'MinG')
set(gca,'Fontsize',14,'LineWidth',1)
ylabel('G max and min')


subplot(2,1,2);
plot(sensfac,maxG-minG,'b*-')
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('k2 parameter scaling')
ylabel('G amplitude (maxG-minG)')
sgtitle('K2 Sensitivity')

%% 
k1 = 0.15;
k2 = 0.2;
k3_def = 0.25;
% initial conditions
y0 = [0; 1; 1];

%parameters
t_end = 5000;
t_transient = 4500;

sensfac = 0.1:0.1:2;

maxG = zeros(length(sensfac),1);
minG = zeros(length(sensfac),1);

for i =1:length(sensfac)
    
    k3 = sensfac(i)*k3_def;

    options = odeset('RelTol',1e-8);
    [t,y] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);

    %determine max and min of G, after transient
    mini = min(find(t > t_transient)); 
    l = length(t);
    maxG(i) = max(y(mini:l,3));
    minG(i) = min(y(mini:l,3));
        
end

%plot numerical solution
figure()
subplot(2,1,1);
plot(sensfac,maxG,'b*-');
hold on;
plot(sensfac,minG,'r*-');
legend('MaxG', 'MinG')
set(gca,'Fontsize',14,'LineWidth',1)
ylabel('G max and min')


subplot(2,1,2);
plot(sensfac,maxG-minG,'b*-')
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('k3 parameter scaling')
ylabel('G amplitude (maxG-minG)')
sgtitle('K3 Sensitivity')
%% 
k1 = 0.15;
%k1 = 0.15;
k2 = 0.2;
k3 = 0.25;
k1_plus = k1 + .05*k1;
k1_minus = k1 - .05*k1;
[tbase,ybase] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tbase > t_transient)); 
l = length(tbase);
amp_base = max(ybase(mini:l,3))- min(ybase(mini:l,3));

[tplus,yplus] = ode45(@(t,y) HPG(t,y, k1_plus,k2,k3),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tplus > t_transient)); 
l = length(tplus);
amp_plus = max(yplus(mini:l,3))- min(yplus(mini:l,3));

[tminus,yminus] = ode45(@(t,y) HPG(t,y, k1_minus,k2,k3),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tminus > t_transient)); 
l = length(tminus);
amp_minus = max(yminus(mini:l,3))- min(yminus(mini:l,3));

dp = k1_plus - k1_minus;
ds = amp_plus - amp_minus;
p = k1;
s = amp_base;
k1_local_sensitivity = (p*ds)/(dp*s);

%% 
%% 
k1 = 0.15;
k2 = 0.2;
k3 = 0.25;
k2_plus = k2 + .05*k2;
k2_minus = k2 - .05*k2;
[tbase,ybase] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tbase > t_transient)); 
l = length(tbase);
amp_base = max(ybase(mini:l,3))- min(ybase(mini:l,3));

[tplus,yplus] = ode45(@(t,y) HPG(t,y, k1,k2_plus,k3),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tplus > t_transient)); 
l = length(tplus);
amp_plus = max(yplus(mini:l,3))- min(yplus(mini:l,3));

[tminus,yminus] = ode45(@(t,y) HPG(t,y, k1,k2_minus,k3),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tminus > t_transient)); 
l = length(tminus);
amp_minus = max(yminus(mini:l,3))- min(yminus(mini:l,3));

dp = k2_plus - k2_minus;
ds = amp_plus - amp_minus;
p = k2;
s = amp_base;
k2_local_sensitivity = (p*ds)/(dp*s);


%% 
k1 = 0.15;
k2 = 0.2;
k3 = 0.25;
k3_plus = k3 + .05*k3;
k3_minus = k3 - .05*k3;
[tbase,ybase] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tbase > t_transient)); 
l = length(tbase);
amp_base = max(ybase(mini:l,3))- min(ybase(mini:l,3));

[tplus,yplus] = ode45(@(t,y) HPG(t,y, k1,k2,k3_plus),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tplus > t_transient)); 
l = length(tplus);
amp_plus = max(yplus(mini:l,3))- min(yplus(mini:l,3));

[tminus,yminus] = ode45(@(t,y) HPG(t,y, k1,k2,k3_minus),[0 t_end],y0,options);
%determine max and min of G, after transient
mini = min(find(tminus > t_transient)); 
l = length(tminus);
amp_minus = max(yminus(mini:l,3))- min(yminus(mini:l,3));

dp = k3_plus - k3_minus;
ds = amp_plus - amp_minus;
p = k3;
s = amp_base;
k3_local_sensitivity = (p*ds)/(dp*s);

%% 
k1_def = 0.15;
k2_def = 0.2;
k3 = 0.25;
% initial conditions
y0 = [0; 1; 1];

%parameters
t_end = 5000;
t_transient = 4500;

sensfac = 0.1:0.1:2;

ampG = zeros(length(sensfac),2);

for i =1:length(sensfac)
    for j=1:length(sensfac)
    
        k1 = sensfac(i)*k1_def;
        k2 = sensfac(j)*k2_def;
    
        options = odeset('RelTol',1e-8);
        [t,y] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);
    
        %determine max and min of G, after transient
        mini = min(find(t > t_transient)); 
        l = length(t);
        ampG(i,j) = max(y(mini:l,3))-min(y(mini:l,3));        
    end
end

figure();
surf(sensfac,sensfac,ampG);
colorbar;
xlabel('K1 Scaling')
ylabel('K2 Scaling')
zlabel('Amplitude (Gmax-Gmin)')
title('K1 and K2');


[C,I] = max(ampG(:));
[I1,I2] = ind2sub(size(ampG),I);
bestk1 = sensfac(I1)*k1_def;
bestk2 = sensfac(I2)*k2_def;

%% 
k1_def = 0.15;
k2 = 0.2;
k3_def = 0.25;
% initial conditions
y0 = [0; 1; 1];

%parameters
t_end = 5000;
t_transient = 4500;

sensfac = 0.1:0.1:2;

ampG = zeros(length(sensfac),2);

for i =1:length(sensfac)
    for j=1:length(sensfac)
    
        k1 = sensfac(i)*k1_def;
        k3 = sensfac(j)*k3_def;
    
        options = odeset('RelTol',1e-8);
        [t,y] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);
    
        %determine max and min of G, after transient
        mini = min(find(t > t_transient)); 
        l = length(t);
        ampG(i,j) = max(y(mini:l,3))-min(y(mini:l,3));        
    end
end

figure();
surf(sensfac,sensfac,ampG);
colorbar;
xlabel('K1 Scaling')
ylabel('K3 Scaling')
zlabel('Amplitude (Gmax-Gmin)')
title('K1 and K3');

[C,I] = max(ampG(:));
[I1,I2] = ind2sub(size(ampG),I);



%% 
k1 = 0.15;
k2_def = 0.2;
k3_def = 0.25;
% initial conditions
y0 = [0; 1; 1];

%parameters
t_end = 5000;
t_transient = 4500;

sensfac = 0.1:0.1:2;

ampG = zeros(length(sensfac),2);

for i =1:length(sensfac)
    for j=1:length(sensfac)
    
        k2 = sensfac(i)*k1_def;
        k3 = sensfac(j)*k2_def;
    
        options = odeset('RelTol',1e-8);
        [t,y] = ode45(@(t,y) HPG(t,y, k1,k2,k3),[0 t_end],y0,options);
    
        %determine max and min of G, after transient
        mini = min(find(t > t_transient)); 
        l = length(t);
        ampG(i,j) = max(y(mini:l,3))-min(y(mini:l,3));        
    end
end

figure();
surf(sensfac,sensfac,ampG);
colorbar;
xlabel('K2 Scaling')
ylabel('K3 Scaling')
zlabel('Amplitude (Gmax-Gmin)')
title('K2 and K3');

[C,I] = max(ampG(:));
[I1,I2] = ind2sub(size(ampG),I);
bestk2 = sensfac(I1)*k2_def;
bestk3 = sensfac(I2)*k3_def;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = HPG(t,y, k1, k2, k3)

n = 9;

dydt = [1/(1+y(3)^n) - k1*y(1); y(1) - k2*y(2); y(2) - k3*y(3)];

end