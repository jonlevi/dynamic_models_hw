%% Part 1
t_end = 200;
options = odeset('RelTol',1e-8);
y0 = [1.8; .1];
beta_values = .1:.4:2;
gamma = 1;
nu = 1;
N=2;
figure();
for j=1:length(beta_values)
    beta = beta_values(j);
    [t,ysolution] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);
    
    %plot solution
    % figure()
    % plot(ysolution(:,1),ysolution(:,2),'--','LineWidth',1)
    hold on;
    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('Susceptible')
    ylabel('Infecteds')
    % title('Beta=1, y0 = [1.8; 0.1]')
    figure();
    subplot(3,1,1);
    plot(t,ysolution(:,1),'b-','LineWidth',2)
    title(strcat('Beta=', num2str(beta)))
    hold on
    ylabel('S')
    set(gca,'Fontsize',14,'LineWidth',1)
    subplot(3,1,2);
    plot(t,ysolution(:,2),'b-','LineWidth',2)
    hold on
    ylabel('I')
    set(gca,'Fontsize',14,'LineWidth',1)
    subplot(3,1,3);
    plot(t,N-ysolution(:,1)-ysolution(:,2),'b-','LineWidth',2)
    hold on
    ylabel('R')
    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('Time (t.u.)')
end
% legend(strcat('\beta=',string(num2cell(beta_values))))


%% Part 2
[sgrid,igrid] = meshgrid(0:0.1:2.2,0:0.1:1);
beta = 1;
gamma = 1;
nu = 1;
N=2;
sdot = -beta*sgrid .* igrid + nu*(N - sgrid - igrid);
idot = beta* sgrid .* igrid - nu*igrid;
figure()
quiver(sgrid, igrid, sdot, idot)
hold on;
plot(2, 0, 'ro')
plot(1, .5, 'ro')
title('Quiver Plot with Fixed Points')
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Susceptible')
ylabel('Infecteds')


figure()
quiver(sgrid, igrid, sdot, idot)
hold on;
plot(2, 0, 'ro')
plot(1, .5, 'ro')
%parameters
t_end = 20;

options = odeset('RelTol',1e-8);

% initial conditions
s_0 = 0:0.1:2;
i_0 = 0:.2:.5;
for idxs=1:length(s_0)
    for idxi=1:length(i_0)
        [t,ysolution] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],[s_0(idxs), i_0(idxi)],options);

        %plot solution
        plot(ysolution(:,1),ysolution(:,2),'r--','LineWidth',1)
        hold on;
        set(gca,'Fontsize',14,'LineWidth',1)
        xlabel('Susceptible')
        ylabel('Infecteds')
        title('Phase Portrait')


    end
end

%% Part 3 
% increase beta
t_end = 20;
options = odeset('RelTol',1e-8);

y0 = [1.8; 0.1];
beta = 1;
gamma = 1;
nu = 1;
N=2;
[t1,ysolution1] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);


beta = 2;
gamma = 1;
nu = 1;
N=2;
[t2,ysolution2] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);
beta = 0.5;
[t3,ysolution3] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);

% plot solution
figure()
plot(ysolution1(:,1),ysolution1(:,2),'r--','LineWidth',1)
hold on;
plot(ysolution2(:,1),ysolution2(:,2),'b--','LineWidth',1)
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Susceptible')
ylabel('Infecteds')
title('Varying Beta, y0 = [1.8; 0.1]')
legend('Beta=1', 'Beta=2')

figure()
plot(ysolution1(:,1),ysolution1(:,2),'r--','LineWidth',1)
hold on;
plot(ysolution3(:,1),ysolution3(:,2),'g--','LineWidth',1)
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Susceptible')
ylabel('Infecteds')
title('Varying Beta, y0 = [1.8; 0.1]')
legend('Beta=1', 'Beta=0.5')



figure()
subplot(3,1,1);
plot(t1,ysolution1(:,1),'r--','LineWidth',2);
hold on;
plot(t2,ysolution2(:,1),'b--','LineWidth',2)
legend('Beta=1', 'Beta=2')
ylabel('S')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(t1,ysolution1(:,2),'r--','LineWidth',2);
hold on;
plot(t2,ysolution2(:,2),'b--','LineWidth',2);
legend('Beta=1', 'Beta=2')
ylabel('I')

set(gca,'Fontsize',14,'LineWidth',1)

subplot(3,1,3);
plot(t1,N-ysolution1(:,1)-ysolution1(:,2),'r--','LineWidth',2)
hold on;
plot(t2,N-ysolution2(:,1)-ysolution2(:,2),'b--','LineWidth',2)

legend('Beta=1', 'Beta=2')


ylabel('R')
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)')
sgtitle('Time Courses for varied Beta')

figure()
subplot(3,1,1);
plot(t1,ysolution1(:,1),'r--','LineWidth',2);
hold on;
plot(t3,ysolution3(:,1),'g--','LineWidth',2)
legend('Beta=1', 'Beta=0.5')
ylabel('S')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(t1,ysolution1(:,2),'r--','LineWidth',2);
hold on;
plot(t3,ysolution3(:,2),'g--','LineWidth',2);
legend('Beta=1', 'Beta=0.5')
ylabel('I')

set(gca,'Fontsize',14,'LineWidth',1)

subplot(3,1,3);
plot(t1,N-ysolution1(:,1)-ysolution1(:,2),'r--','LineWidth',2)
hold on;
plot(t3,N-ysolution3(:,1)-ysolution3(:,2),'g--','LineWidth',2)

legend('Beta=1', 'Beta=0.5')


ylabel('R')
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)')
sgtitle('Time Courses for varied Beta')


%% N = 1

beta = 2;
gamma = 1;
nu = 1;
N=2;
[t1,ysolution1] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);
beta = 2;
gamma = 1;
nu = 1;
N=1;
[t2,ysolution2] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);

% plot solution
figure()
plot(ysolution1(:,1),ysolution1(:,2),'r--','LineWidth',1)
hold on;
plot(ysolution2(:,1),ysolution2(:,2),'b--','LineWidth',1)
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Susceptible')
ylabel('Infecteds')
title('Vaccinaton: Varying N, y0 = [1.8; 0.1]')
legend('N=2', 'N=1')



figure()
subplot(3,1,1);
plot(t1,ysolution1(:,1),'r--','LineWidth',2);
hold on;
plot(t2,ysolution2(:,1),'b--','LineWidth',2)
legend('N=2', 'N=1')
ylabel('S')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(t1,ysolution1(:,2),'r--','LineWidth',2);
hold on;
plot(t2,ysolution2(:,2),'b--','LineWidth',2);
legend('N=2', 'N=1')
ylabel('I')

set(gca,'Fontsize',14,'LineWidth',1)

subplot(3,1,3);
plot(t1,N-ysolution1(:,1)-ysolution1(:,2),'r--','LineWidth',2)
hold on;
plot(t2,N-ysolution2(:,1)-ysolution2(:,2),'b--','LineWidth',2)

legend('N=2', 'N=1')


ylabel('R')
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)')
sgtitle('Time Courses for varied N')
%% %% nu = 2

beta = 2;
gamma = 1;
nu = 1;
N=2;
[t1,ysolution1] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);
beta = 2;
gamma = 1;
nu = 2;
N=2;
[t2,ysolution2] = ode45(@(t,y) [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)],[0 t_end],y0,options);

% plot solution
figure()
plot(ysolution1(:,1),ysolution1(:,2),'r--','LineWidth',1)
hold on;
plot(ysolution2(:,1),ysolution2(:,2),'b--','LineWidth',1)
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Susceptible')
ylabel('Infecteds')
title('Testing/Stay at Home: Varying \nu, y0 = [1.8; 0.1], \beta=2')
legend('\nu=1', '\nu=2')



figure()
subplot(3,1,1);
plot(t1,ysolution1(:,1),'r--','LineWidth',2);
hold on;
plot(t2,ysolution2(:,1),'b--','LineWidth',2)
legend('\nu=1', '\nu=2')
ylabel('S')
set(gca,'Fontsize',14,'LineWidth',1)


subplot(3,1,2);
plot(t1,ysolution1(:,2),'r--','LineWidth',2);
hold on;
plot(t2,ysolution2(:,2),'b--','LineWidth',2);
legend('\nu=1', '\nu=2')
ylabel('I')

set(gca,'Fontsize',14,'LineWidth',1)

subplot(3,1,3);
plot(t1,N-ysolution1(:,1)-ysolution1(:,2),'r--','LineWidth',2)
hold on;
plot(t2,N-ysolution2(:,1)-ysolution2(:,2),'b--','LineWidth',2)

legend('\nu=1', '\nu=2')


ylabel('R')
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)')
sgtitle('Time Courses for varied N')

