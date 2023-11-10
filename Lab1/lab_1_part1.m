% initial condition
x0 = 10;

%parameters
t_end = 5;

options = odeset('RelTol',1e-8);
[t,x] = ode45(@exp_decay,[0 t_end],x0,options);


dt = t_end/(length(t)-1);
t_even = 0:dt:t_end;

%plot numerical solution
hold on
plot(t,x,'b.-')
plot(t_even, x, 'r.-')
legend('With t from simulation','With t evenly spaced')
xlabel('t'), ylabel('x')
title('Solving dx/dt = -kx')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dxdt = exp_decay(t,x)

k = 2;

dxdt = -k*x;

end
