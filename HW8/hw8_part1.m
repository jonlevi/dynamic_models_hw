
noise_levels = 0:500:10000;
% noise_levels = logspace(2.8,4, 10);
% noise_levels = [5000];
Nsteps = 20000;
mean_values = zeros(length(noise_levels),1);
variance_values = zeros(length(noise_levels),1);
distributions = zeros(length(noise_levels), 20);
for n=1:length(noise_levels)

    
    a = 0; 
    b = 0.8;
    c = 3.0;
    
    
    y1_0 = -1;
    y2_0 = -1;
    
    
    dt=0.01;
    Dnz=noise_levels(n);
    
    t = zeros(1, Nsteps);
    
    y1 = y1_0*ones(1,Nsteps);
    y2 = y2_0*ones(1,Nsteps);
    sig = sqrt(2*Dnz);
    fac = sig*sqrt(dt);
    S=0;
    for i=2:Nsteps
        t(i) = (i-1)*dt;
    
        dy1 = c*(y1(i-1)-1/3*(y1(i-1)).^3-y2(i-1)+S) + fac*randn;
        dy2 = (y1(i-1)+a-b*y2(i-1))/c;
        y1(i) = y1(i-1)+dy1*dt;
        y2(i) = y2(i-1)+dy2*dt;
    end
    
    % %parameters
    % t_end = 50;
    % 
    % options = odeset('RelTol',1e-8,'Maxstep',0.1);
    % [t,y] = ode15s(@(t,y) FHNeqns(t,y,params),[0 t_end],y0,options);
    % %plot numerical solution


    [Voltage_pks,locs] = findpeaks(y1,'MinPeakDistance',500,'MinPeakProminence',1);
    % 
    % figure(2);
    % figure(n+1)
    % plot(t,y1,'LineWidth',1)
    % title(num2str(Dnz))
    % hold on
    % plot(t(locs), Voltage_pks, 'ko')
    % ylabel('x')
    % set(gca,'Fontsize',14,'LineWidth',1);
    % subplot(212);
    % plot(t,y2,'LineWidth',1)
    % hold on
    % set(gca,'Fontsize',14,'LineWidth',1)
    % xlabel('Time (t.u.)'), ylabel('y')

    % figure(2);
    % subplot(6,1,n);
    % hold on;
    % histogram(diff(t(locs)),8.5:.5:10,'Normalization','probability');
    m = mean(diff(t(locs)));
    v = var(diff(t(locs)));
    isi = diff(t(locs));
    distributions(n,:) = isi(1:20);
    mean_values(n) = m;
    variance_values(n) = v;
    % plot([m m] ,[0 1], 'r--')

    

    % 
    % %plot numerical solution
    % figure(2)
    % plot(y1,y2,'LineWidth',2)
    % hold on
    % set(gca,'Fontsize',14,'LineWidth',1)
    % xlabel('x'), ylabel('y')

    % figure(2);
    % subplot(5,2,2*n-1);
    % histHandle = histogram(y1(t>20)+1.2, -.2:.05:.2, 'Normalization','probability');
    % histHandle.BinEdges = histHandle.BinEdges + histHandle.BinWidth/2;
    % title(strcat('D=', num2str(Dnz)))
    % ylabel('x')
    % subplot(5,2,2*n);
    % histHandle = histogram(y2(t>20)+.62, -.06:.01:.06,'Normalization','probability');
    % histHandle.BinEdges = histHandle.BinEdges + histHandle.BinWidth/2;
    % title(strcat('D=', num2str(Dnz)))
    % ylabel('y')
   
end
figure(1);


boxplot(distributions', noise_levels,'PlotStyle','compact', 'LabelOrientation','horizontal');
xlabel('Noise Intensity (D)', 'FontSize',14);
ylabel('Interspike Interval', 'FontSize',14);
xtickangle(45)

figure(2);
x = noise_levels;
y = mean_values;
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(x, y, 'b.', 'MarkerSize', 15);
hold on;
plot(xFit, yFit, 'r--', 'LineWidth', 2); 
xlabel('Noise Intensity (D)')
ylabel('Mean ISI')

figure(3);
x = noise_levels;
y = variance_values;
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(x, y, 'b.', 'MarkerSize', 15);
hold on;
plot(xFit, yFit, 'r--', 'LineWidth', 2); 
xlabel('Noise Intensity (D)')
ylabel('Variance if ISI')
% bar(noise_levels, mean_values)
% xlim([min(log10(noise_levels)), max(log10(noise_levels))])
% set(gca,'XScale','log');
% figure(2);
% subplot(5,2,2*n-1);
% xlabel('x(t)-x*','FontSize',16);
% subplot(5,2,2*n);
% xlabel('y(t)-y*', 'FontSize',16);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
noise_levels = 0:.1:.5;
% noise_levels = logspace(2.8,4, 10);
% noise_levels = [5000];
Nsteps = 15000;
mean_values = zeros(length(noise_levels),1);
variance_values = zeros(length(noise_levels),1);
distributions = zeros(length(noise_levels), 20);
for n=1:length(noise_levels)

    
    a = .7; 
    b = 0.8;
    c = 3.0;
    
    
    y1_0 = -1;
    y2_0 = -1;
    
    
    dt=0.01;
    Dnz=noise_levels(n);

    cn = dsp.ColoredNoise('Color','pink','SamplesPerFrame',Nsteps);
    noise = cn();
    noise = noise/std(noise);
    t = zeros(1, Nsteps);
    
    y1 = y1_0*ones(1,Nsteps);
    y2 = y2_0*ones(1,Nsteps);
    sig = sqrt(2*Dnz);
    fac = sig*sqrt(dt);
    S=0;
    for i=2:Nsteps
        t(i) = (i-1)*dt;
    
        dy1 = c*(y1(i-1)-1/3*(y1(i-1)).^3-y2(i-1)+S);
        dy2 = (y1(i-1)+a-b*y2(i-1))/c + fac*noise(i-1);
        y1(i) = y1(i-1)+dy1*dt;
        y2(i) = y2(i-1)+dy2*dt;
    end
    figure(1);
    subplot(2,1,1);
    plot(t,y1,'LineWidth',1)
    hold on
    % plot(t(locs), Voltage_pks, 'ko')
    ylabel('x')
    set(gca,'Fontsize',14,'LineWidth',1);
    subplot(212);
    plot(t,y2,'LineWidth',1)
    hold on
    set(gca,'Fontsize',14,'LineWidth',1)
    xlabel('Time (t.u.)'), ylabel('y')


end
subplot(2,1,1);
legend(strcat('D=',string(num2cell(noise_levels))))

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


function dydt = FHNeqns(t,y, params)


dydt = [params.c*(y(1)-1/3*y(1)^3-y(2)); (y(1)+params.a-params.b*y(2))/params.c];

end