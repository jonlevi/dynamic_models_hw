
global D;
global L;
global k;
global a;


D=10;
k=1;
L=10;
a=1;
t_end = 10;

%set up solution mesh
m = 0;
Nx = 200;
Nt = 300;
x = linspace(0,L,Nx);
t = linspace(0,t_end,Nt);



D=8;
%solve the partial differential equations
sol = pdepe(m,@rf_eqn,@rf_initial,@rf_bc,x,t);


C = sol(:,:);

% plot results
figure(1)
subplot(3,1,1);
surf(x,t,C, 'EdgeColor','none')    
xlabel('Position (\mum)')
ylabel('Time (s)')
zlabel('C(x,t)')
shading interp
title(strcat('D= ', num2str(D)))
set(gca,'fontsize',14)
colorbar;
clim([0, .3]);
hold on;
view(2);
% 
figure(2);
plot(x, C(end,:))
xlabel('x')
ylabel('C(x,10)')

D=10;
sol = pdepe(m,@rf_eqn,@rf_initial,@rf_bc,x,t);

C = sol(:,:);

% % plot results
figure(1)
subplot(3,1,2);
surf(x,t,C, 'EdgeColor','none')    
xlabel('Position (\mum)')
ylabel('Time (s)')
zlabel('C(x,t)')
shading interp
title(strcat('D= ', num2str(D)))
set(gca,'fontsize',14)
view(2)
colorbar;
clim([0, .3]);
% clim([-60, 35]);
hold on;

figure(2);
hold on;
plot(x, C(end,:))
xlabel('x')
ylabel('C(x,10)')

D=12;
sol = pdepe(m,@rf_eqn,@rf_initial,@rf_bc,x,t);

C = sol(:,:);
% 
% % % plot results
figure(1)
subplot(3,1,3);
surf(x,t,C, 'EdgeColor','none')    
xlabel('Position (\mum)')
ylabel('Time (s)')
zlabel('C(x,t)')
shading interp
title(strcat('D= ', num2str(D)))
set(gca,'fontsize',14)
view(2)
colorbar;
clim([0, .3]);
hold on;
% 
figure(2);
plot(x, C(end,:))
xlabel('x')
ylabel('C(x,10)')
legend('D=8', 'D=10', 'D=12')
%% 
D=10;
L_values = 5:1:25;
left_values = zeros(1,length(L_values));
right_values = zeros(1,length(L_values));
for i=1:length(L_values)
    L = L_values(i);
    x = linspace(0,L,Nx);
    sol = pdepe(m,@rf_eqn,@rf_initial,@rf_bc,x,t);
    C = sol(:,:);

    % figure(2);
    % plot(x, C(end,:))
    % hold on;
    % xlabel('x')
    % ylabel('C(x,10)')

    left_values(i) = C(end,1);
    right_values(i) = C(end, end);

end
% figure(2);
% legend(strcat('L',string(num2cell(L_values))))
figure(3);
subplot(2,1,1);
plot(L_values, left_values);
ylabel('C(0,10)')
subplot(2,1,2);
plot(L_values, right_values);
xlabel('L')
ylabel('C(L,10)')


figure(4);
plot(L_values, left_values);
ylabel('C(0,10)')
xlabel('L');
hold on;
plot([min(L_values) max(L_values)], [a/sqrt(k*D) a/sqrt(k*D)], 'r--')

%% 


%set equations
function [c,f,s] = rf_eqn(x,t,u,DuDx)

global D;
global k;

c = 1;
f = D*DuDx;
s = -1*k*u; 

end


%set initial condition
function value =rf_initial(x)
value=0;
end


%set boundary conditions
% form is p + q(D dudx)=0
function [pl,ql,pr,qr] = rf_bc(xl,ul,xr,ur,t)

global D;
global a;

% Left boundary x=0; dv/dx=-a
ql = 1;
pl=a;

% Right boundary x=L; dv/dx=0
pr=0;
qr=1;


end