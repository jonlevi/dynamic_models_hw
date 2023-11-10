% initial conditions
y0 = [1.8; 0.1];

%parameters
t_end = 20;
global N;
N = 2;

options = odeset('RelTol',1e-8);
[t,y] = ode45(@SIeqns,[0 t_end],y0,options);

%plot solution
figure(1)
plot(y(:,1),y(:,2),'b-','LineWidth',2)
hold on
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Susceptible')
ylabel('Infecteds')

figure(2)
subplot(311), plot(t,y(:,1),'b-','LineWidth',2)
hold on
ylabel('S')
set(gca,'Fontsize',14,'LineWidth',1)
subplot(312), plot(t,y(:,2),'b-','LineWidth',2)
hold on
ylabel('I')
set(gca,'Fontsize',14,'LineWidth',1)
subplot(313), plot(t,N-y(:,1)-y(:,2),'b-','LineWidth',2)
hold on
ylabel('R')
set(gca,'Fontsize',14,'LineWidth',1)
xlabel('Time (t.u.)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = SIeqns(t,y)

global N;
beta = 1;
gamma = 1;
nu = 1;

dydt = [-beta*y(1)*y(2) + gamma*(N-y(1)-y(2)); beta*y(1)*y(2) - nu*y(2)];

end
