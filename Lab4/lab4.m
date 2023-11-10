lag = 0.9;

t_end = 120;

opts = ddeset('RelTol',1e-8,'AbsTol',1e-8);
initial_x = .1:.5:1;
initial_y = .1:.5:1;
initial_s = .1:.5:1;
for x_ind = 1:length(initial_x)
    % for y_ind = 1:length(initial_x)
    %     for s_ind = 1:length(initial_x)
            sol = dde23(@p53de,lag,@(t) ...
                [initial_x(x_ind) initial_x(y_ind) initial_x(s_ind)],[0, t_end],opts);

            figure(1);
            hold on;
            subplot(3,1,1);
            plot(sol.x,sol.y(1,:),'--','LineWidth',1)
            hold on;
            ylabel('X')
            set(gca,'Fontsize',14,'LineWidth',1)
            
            subplot(3,1,2);
            plot(sol.x,sol.y(2,:),'--','LineWidth',1);
            hold on;
            ylabel('Y')
            set(gca,'Fontsize',14,'LineWidth',1)
            subplot(3,1,3);
            plot(sol.x,sol.y(3,:),'--','LineWidth',1);
            hold on;
            hold on, xlabel('Time (hours)'), ylabel('S');
            set(gca,'Fontsize',14,'LineWidth',1)
            
            figure(2);
            plot(sol.y(1,:),sol.y(2,:),'--','LineWidth',1)
            xlabel('p53');
            ylabel('Mdm2');
            set(gca,'Fontsize',14,'LineWidth',1)
            hold on;
            
            figure(3);
            plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'--','LineWidth',1)
            xlabel('p53');
            ylabel('Mdm2');
            zlabel('S');
            set(gca,'Fontsize',14,'LineWidth',1)
            hold on;
           
    %     end
    % end
end
figure()
plot(sol.x,sol.y(1,:),'r--','LineWidth',1);
hold on;
plot(sol.x,sol.y(2,:),'b--','LineWidth',1);
legend('X', 'Y')
title('Model with S')
xlabel('Time (hours)')

%% 
t_end = 500;
lag = 0.9;
sol = dde23(@p53deChangedN,lag,@(t) ...
                [.1,.1,.5],[0, t_end],opts);
figure(1);
title('n=2')
hold on;
subplot(3,1,1);
plot(sol.x,sol.y(1,:),'r-','LineWidth',1)
hold on;
ylabel('X')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(sol.x,sol.y(2,:),'r-','LineWidth',1);
hold on;
ylabel('Y')
set(gca,'Fontsize',14,'LineWidth',1)
subplot(3,1,3);
plot(sol.x,sol.y(3,:),'r-','LineWidth',1);
hold on;
hold on, xlabel('Time (hours)'), ylabel('S');
set(gca,'Fontsize',14,'LineWidth',1)

figure(2);
title('n=2')
plot(sol.y(1,:),sol.y(2,:),'r-','LineWidth',1)
xlabel('p53');
ylabel('Mdm2');
set(gca,'Fontsize',14,'LineWidth',1)
hold on;

figure(3);
title('n=2')
plot3(sol.y(1,:),sol.y(2,:),sol.y(3,:),'r-','LineWidth',1)
xlabel('p53');
ylabel('Mdm2');
zlabel('S');
set(gca,'Fontsize',14,'LineWidth',1)
hold on;
%% 
t_end=400;
sol = dde23(@p53deChangedalpha,lag,@(t) ...
                [.1,.1,.5],[0, t_end],opts);
figure(1);
sgtitle('\alpha_{xy} = 0.5')
hold on;
subplot(3,1,1);
plot(sol.x,sol.y(1,:),'r-','LineWidth',1)
hold on;
ylabel('X')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(sol.x,sol.y(2,:),'r-','LineWidth',1);
hold on;
ylabel('Y')
set(gca,'Fontsize',14,'LineWidth',1)
subplot(3,1,3);
plot(sol.x,sol.y(3,:),'r-','LineWidth',1);
hold on;
hold on, xlabel('Time (hours)'), ylabel('S');
set(gca,'Fontsize',14,'LineWidth',1)
%% %% 
t_end=400;
sol = dde23(@p53de,0.5,@(t) ...
                [.1,.1,.5],[0, t_end],opts);
figure(1);
sgtitle('\tau = 0.5')
hold on;
subplot(3,1,1);
plot(sol.x,sol.y(1,:),'r-','LineWidth',1)
hold on;
ylabel('X')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(sol.x,sol.y(2,:),'r-','LineWidth',1);
hold on;
ylabel('Y')
set(gca,'Fontsize',14,'LineWidth',1)
subplot(3,1,3);
plot(sol.x,sol.y(3,:),'r-','LineWidth',1);
hold on;
hold on, xlabel('Time (hours)'), ylabel('S');
set(gca,'Fontsize',14,'LineWidth',1)
%% 
t_end=50;
lag = 6.3;
sol = dde23(@p53AlternativeModel,lag,@(t) ...
                [.15,.2],[0, t_end],opts);
sol2 = dde23(@p53AlternativeModel,lag,@(t) ...
                [.1,.34],[0, t_end],opts);
figure(1);

sgtitle('Alternative Model Without s')
hold on;
subplot(3,1,1);

plot(sol.x,sol.y(1,:),'r-','LineWidth',1)
ylim([0.05 .25])
hold on;
plot(sol2.x,sol2.y(1,:),'b-','LineWidth',1)
ylabel('X')
legend('x_0, y_0 = (.15,.20)' ,'x_0, y_0 = (.10,.34)')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(sol.x,sol.y(2,:),'r-','LineWidth',1);
hold on;
plot(sol2.x,sol2.y(2,:),'b-','LineWidth',1);
ylim([0.05 .35])
hold on;
ylabel('Y')
set(gca,'Fontsize',14,'LineWidth',1)

figure()

plot(sol.x,sol.y(1,:),'r--','LineWidth',1);
hold on;
plot(sol.x,sol.y(2,:),'b--','LineWidth',1);
ylim([0.1 .17])
legend('X', 'Y')
title('Model without S')
xlabel('Time (hours)')

%% %% 
t_end=50;
lag = 6.3;
sol = dde23(@p53AlternativeModel,lag,@(t) ...
                [.2,.2],[0, t_end],opts);
sol2 = dde23(@p53AlternativeModelDoubleBeta,lag,@(t) ...
                [.2,.2],[0, t_end],opts);
figure(1);

sgtitle('Alternative Model')
hold on;
subplot(3,1,1);

plot(sol.x,sol.y(1,:),'r-','LineWidth',1)
hold on;
plot(sol2.x,sol2.y(1,:),'b-','LineWidth',1)
ylim([0, .3])
ylabel('X')
legend('\beta_x=2.3' ,'\beta_x=4.6')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(sol.x,sol.y(2,:),'r-','LineWidth',1);
ylim([0, .3])
hold on;
plot(sol2.x,sol2.y(2,:),'b-','LineWidth',1);
hold on;
ylabel('Y')
set(gca,'Fontsize',14,'LineWidth',1)

% figure()
% 
% plot(sol.x,sol.y(1,:),'r--','LineWidth',1);
% hold on;
% plot(sol.x,sol.y(2,:),'b--','LineWidth',1);
% ylim([0.1 .17])
% legend('X', 'Y')
% title('Model without S, with double beta')
% xlabel('Time (hours)')
%% %% %% 
t_end=200;
lag = 0.9;
sol = dde23(@p53de,lag,@(t) ...
                [.1,.1, .1],[0, t_end],opts);
sol2 = dde23(@p53ChangedBeta,lag,@(t) ...
                [.1,.1, .1],[0, t_end],opts);
figure(1);

sgtitle('Original Model (with S)')
hold on;
subplot(3,1,1);

plot(sol.x,sol.y(1,:),'r-','LineWidth',1)
hold on;
plot(sol2.x,sol2.y(1,:),'b-','LineWidth',1)
ylabel('X')
legend('\beta_x=0.9' ,'\beta_x=1.8')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(sol.x,sol.y(2,:),'r-','LineWidth',1);
hold on;
plot(sol2.x,sol2.y(2,:),'b-','LineWidth',1);
hold on;
ylabel('Y')
set(gca,'Fontsize',14,'LineWidth',1)

%% 


% --------------------------------------------------------------------------

function dydt = p53de(t,y,Z)
% Differential equations function 
betay = 1.2;
alphay = 0.8;
betax = 0.9;
alphaxy = 1.4;
betas = 0.9;
alphas = 2.7;
n = 4;

ylag = Z(:,1);
dydt = [ betax*((y(3))^n)/(1+(y(3))^n) - alphaxy*y(2)*y(1);
         betay*ylag(1) - alphay*y(2);
         betas - alphas*y(2)*y(3)];

end
function dydt = p53ChangedBeta(t,y,Z)
% Differential equations function 
betay = 1.2;
alphay = 0.8;
betax = 2*0.9;
alphaxy = 1.4;
betas = 0.9;
alphas = 2.7;
n = 4;

ylag = Z(:,1);
dydt = [ betax*((y(3))^n)/(1+(y(3))^n) - alphaxy*y(2)*y(1);
         betay*ylag(1) - alphay*y(2);
         betas - alphas*y(2)*y(3)];

end

function dydt = p53deChangedN(t,y,Z)
% Differential equations function 
betay = 1.2;
alphay = 0.8;
betax = 0.9;
alphaxy = 2;
betas = 0.9;
alphas = 2.7;
n = 2;

ylag = Z(:,1);
dydt = [ betax*((y(3))^n)/(1+(y(3))^n) - alphaxy*y(2)*y(1);
         betay*ylag(1) - alphay*y(2);
         betas - alphas*y(2)*y(3)];

end

function dydt = p53deChangedalpha(t,y,Z)
% Differential equations function 
betay = 1.2;
alphay = 0.8;
betax = 0.9;
alphaxy = .5;
betas = 0.9;
alphas = 2.7;
n = 4;

ylag = Z(:,1);
dydt = [ betax*((y(3))^n)/(1+(y(3))^n) - alphaxy*y(2)*y(1);
         betay*ylag(1) - alphay*y(2);
         betas - alphas*y(2)*y(3)];

end

function dydt = p53AlternativeModel(t,y,Z)
betay = 24;
alphay = 24;
betax = 2.3;
alphaxy = 120;
ylag = Z(:,1);
dydt = [ betax - alphaxy*y(2)*y(1);
         betay*ylag(1) - alphay*y(2);];

end
function dydt = p53AlternativeModelDoubleBeta(t,y,Z)
betay = 24;
alphay = 24;
betax = 2.3*2;
alphaxy = 120;
ylag = Z(:,1);
dydt = [ betax - alphaxy*y(2)*y(1);
         betay*ylag(1) - alphay*y(2);];

end