load tumor_growth_data.txt
xdata = tumor_growth_data(:,1);
ydata = tumor_growth_data(:,2);

% fitting and plotting
fun = @(x) x(1)./(x(2)+exp(-x(3)*xdata))-ydata;
x0 = [.4,.5, .4];   %initital guess
options = optimoptions('lsqnonlin','Display','iter');
options.Algorithm = 'levenberg-marquardt';
xf = lsqnonlin(fun,x0,[],[],options)

subplot(3,1,1)
plot(xdata,ydata,'ko')
hold on
xvals = linspace(min(xdata),max(xdata),100);
plot(xvals,xf(1)./(xf(2)+exp(-xf(3).*xvals)),'b-','LineWidth',2)
set(gca,'Fontsize',12,'LineWidth',1)
title('Data and fit')
xlabel('Time')
ylabel('Growth')

% residuals
r = ydata - (xf(1)./(xf(2)+exp(-xf(3).*xdata)));
subplot(3,1,2)
plot(xdata,r,'o')
set(gca,'Fontsize',12,'LineWidth',1)
title('Residuals')
xlabel('Time')
ylabel('Residuals')

subplot(3,2,5)
hist(r)
set(gca,'Fontsize',12,'LineWidth',1)
title('Residuals histogram')
xlabel('Residuals')
ylabel('Number')

% qqplot
subplot(3,2,6)
qqplot(r)
set(gca,'Fontsize',12,'LineWidth',1)
title('QQ plot'), box on
xlabel('Theoretical quantiles')
ylabel('Observed quantiles')


k1 = xf(3);
k2 = xf(2)*xf(3)/xf(1);
options = odeset('RelTol',1e-8);
t_end = 400;

[t,ysolution] = ode45(@(t,x) k1*x-k2*x.^2,[0 t_end],2,options);
figure()
plot(t, ysolution, 'LineWidth' ,2);
hold on;
plot(xvals,xf(1)./(xf(2)+exp(-xf(3).*xvals)),'k--','LineWidth',.5)
plot(xdata,ydata,'ko');
legend('Model Integration', 'Least Squares Logistic Fit', 'Actual Data')
xlabel('time')
ylabel('x(t)')
set(gca,'Fontsize',14)

