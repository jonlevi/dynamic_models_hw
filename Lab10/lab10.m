%stochastic implementation (Gillespie SSA) of constitutive gene expression
%and mRNA degradation

gr_values = .5:.5:5.5;
kr_values = 5:1:15;
mean_values = zeros(length(gr_values),1);
variance_values = zeros(length(gr_values),1);
% distributions = zeros(length(gr_values), 5000);

% figure(2);
% layout = tiledlayout('flow');


burstFactors = 1:5:20;
nLoop=length(burstFactors);
for index=1:nLoop
    burstFactor=burstFactors(index);

    
    %set kinetic parameter values
    kr=10/burstFactor;
    gr=1;
    
    %set initial condition
    X=0;
    j=1;
    t=0;
    
    
    while t < 700
    
        %calculate propensities
        a1=kr;
        a2=gr*X(j);
        asum=a1+a2;
    
        %update counter
        j=j+1;
    
        %update time
        t(j)=t(j-1)-log(rand(1))/asum;
        
        %state transition
        mu=rand(1);
        if mu < a1/asum
            X(j)=X(j-1)+burstFactor;
        else
            X(j)=max(X(j-1)-1,0);
        end
    
    end
    
    % figure(1)
    % stairs(t, X, 'b', 'linewidth', 2)
    % ylabel('Number of mRNA molecules')
    % xlabel('Time (t.u.)')
    % title('Burst=5')
    % set(gca,'fontsize',14)
    % hold on
    % 
    non_transient = X(end-4999:end);
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
    figure(2);
    hold on;
    pspectrum(non_transient);
    disp(var(non_transient))

end
legend(strcat('BurstFactor=',string(num2cell(burstFactors))));
title('Power Spectrum of non-transient M(t)')
xlabel('Frequency')
ylabel('Power (a.u)')
% ylabel(layout, 'Number of events')
% 
% xlabel(layout, 'Number of mRNA molecules')

% figure();
% subplot(1,2,1);
% x = 1./gr_values;
% y = mean_values;
% coefficients = polyfit(x, y, 1);
% xFit = linspace(min(x), max(x), 1000);
% yFit = polyval(coefficients , xFit);
% plot(x, y, 'b.', 'MarkerSize', 15);
% hold on;
% plot(xFit, yFit, 'r--', 'LineWidth', 2); 
% xlabel('1/dr')
% ylabel('Mean value of M')
% set(gca,'fontsize',15)
% subplot(1,2,2);
% x = -1*log10(gr_values);
% y = log10(variance_values);
% coefficients = polyfit(x, y, 1);
% xFit = linspace(min(x), max(x), 1000);
% yFit = polyval(coefficients , xFit);
% plot(x, y, 'b.', 'MarkerSize', 15);
% hold on;
% plot(xFit, yFit, 'r--', 'LineWidth', 2); xlabel('-log(dr)')
% ylabel('log(COV)')
% set(gca,'fontsize',15)