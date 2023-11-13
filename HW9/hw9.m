%stochastic implementation (Gillespie SSA) of constitutive gene expression
%and mRNA degradation

% gr_values = .5:.5:5.5;
kr_values = 10:10:100;
mean_values = zeros(length(kr_values),2);
variance_values = zeros(length(kr_values),2);
mean_dt_values = zeros(length(kr_values),1);
% distributions = zeros(length(gr_values), 5000);

% figure(2);
% layout = tiledlayout('flow');


% burstFactors = 1:5:20;
% nLoop=length(burstFactors);
nLoop=length(kr_values);
for index=1:nLoop
    % burstFactor=burstFactors(index);

    
    %set kinetic parameter values
    kr = kr_values(index);
    kp = 6;
    dr = 1;
    dp = 1;
    
    %set initial condition
    M=0;
    P=0;
    j=1;
    t=0;
    
    
    while max(t) < 700
    
        %calculate propensities
        a1=kr;
        a2=dr*M(j);
        a3=kp*M(j);
        a4=dr*P(j);
        asum=a1+a2+a3+a4;
    
        %update counter
        j=j+1;
    
        %update time

        t(j)=t(j-1)-log(rand(1))/asum;
        
        %state transition
        mu=rand(1);
        if mu < a1/asum
            M(j)=M(j-1)+1;
            P(j) = P(j-1);
        elseif mu < (a1/asum + a2/asum)
            M(j)=max(M(j-1)-1,0);
            P(j) = P(j-1);
        elseif mu < (a1/asum + a2/asum + a3/asum)
            M(j) = M(j-1);
            P(j) = P(j-1)+1;
        else
            M(j) = M(j-1);
            P(j) = max(P(j-1)-1,0);
        end
    
    end
    mean_dt_values(index) = mean(diff(t));
    
    % figure(1)
    % subplot(2,1,1);
    % stairs(t, M, 'b', 'linewidth', 2)
    % ylabel('N_m')
    % set(gca,'fontsize',14)
    % xlim([0 700])
    % hold on;
    % subplot(2,1,2);
    % stairs(t, P, 'b', 'linewidth', 2)
    % xlim([0 700])
    % ylabel('N_p')
    % xlabel('Time (t.u.)')
    % set(gca,'fontsize',14)

    % 
    % figure(2);
    % subplot(2,1,1);
    non_transient_M = M(1000:end);
    mean_M = mean(non_transient_M);
    cv_M = std(non_transient_M)/mean_M;
    mean_values(index,1) = mean_M;
    variance_values(index,1) = cv_M;
    % histogram(non_transient_M,'FaceAlpha',.4,'FaceColor','b');
    % hold on;
    % plot([mean(non_transient_M) mean(non_transient_M)], [0 7000], 'r--')
    % xlabel('M value')
    % ylabel('Count')
    % 
    % subplot(2,1,2);
    non_transient_P = P(1000:end);
    mean_P = mean(non_transient_P);
    cv_P = std(non_transient_P)/mean_P;
    mean_values(index,2) = mean_P;
    variance_values(index,2) = cv_P;

    % histogram(non_transient_P,'FaceAlpha',.4,'FaceColor','b');
    % hold on;
    % xlabel('P Value')
    % ylabel('Count')
    % plot([mean(non_transient_P) mean(non_transient_P)], [0 2500], 'r--')
    % mean_values(index) = mean(non_transient);
    % cv = std(non_transient)/mean(non_transient);
    % variance_values(index)=cv;
    % 
    % 
    % % nexttile;
    % figure(2);
    % histogram(non_transient,'FaceAlpha',.4,'FaceColor','b')
    % title('Burst=5');
    % hold on;
    % plot([mean(non_transient) mean(non_transient)], [0 900], 'r--')
    % set(gca,'fontsize',14)
    % hold on
    % figure(2);
    % hold on;
    % pspectrum(non_transient);
    % disp(var(non_transient))

end
% legend(strcat('BurstFactor=',string(num2cell(burstFactors))));
% title('Power Spectrum of non-transient M(t)')
% xlabel('Frequency')
% ylabel('Power (a.u)')
% ylabel(layout, 'Number of events')
% 
% xlabel(layout, 'Number of mRNA molecules')

figure();
subplot(2,2,1);
x = kr_values;
y = mean_values(:,1);
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(x, y, 'b.', 'MarkerSize', 15);
hold on;
plot(xFit, yFit, 'r--', 'LineWidth', 2); 
xlabel('kr')
ylabel('Mean value of M')
set(gca,'fontsize',15)
subplot(2,2,2);
x = kr_values;
y = mean_values(:,2);
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(x, y, 'b.', 'MarkerSize', 15);
hold on;
plot(xFit, yFit, 'r--', 'LineWidth', 2);
xlabel('kr')
ylabel('Mean value of P')
set(gca,'fontsize',15)
subplot(2,2,3);
x = log10(kr_values);
y = log10(variance_values(:,1));
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(x, y, 'b.', 'MarkerSize', 15);
hold on;
plot(xFit, yFit, 'r--', 'LineWidth', 2);
xlabel('log(kr)')
ylabel('log(cov) of M')
set(gca,'fontsize',15)


subplot(2,2,4);
x = log10(kr_values);
y = log10(variance_values(:,2));
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(x, y, 'b.', 'MarkerSize', 15);
hold on;
plot(xFit, yFit, 'r--', 'LineWidth', 2);
xlabel('log(kr)')
ylabel('log(cov) of P')
set(gca,'fontsize',15)