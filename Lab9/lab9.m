clear all
n_totals = linspace(1000, 100000, 100);    %number of crossbridges
forces = [400, 40, 4, .4, .04];
% alphas=[10 28];   % /s  probability per unit time for attachment
% betas = [100 180];
colors = ["red" ,"magenta" , "black" ,"blue", "cyan", "#ffa500"];
%  % ,"magenta", "cyan", "#808080", "#8B8000", "#ffa500"];
beta=126;    % /s  probability per unit time for detachment
alpha=14;
% N_a = length(alphas);
% N_b = length(betas);
% N_ex = 5;
% idx=0;
% index = reshape(1:N_b*N_a*N_ex, N_b*N_a, N_ex).';
c_index=0;
cov_values = zeros(length(n_totals),1);

for n_total=1:length(n_totals)
    c_index = c_index+1;
        % alpha = alphas(alp);
        % beta = betas(bet);
    % C = colors(c_index);
    n0 = floor(n_totals(n_total));
    dt=0.01/(alpha+beta);       % s  duration of time step
    tend=0.25;  % simulation time 
    klokmax=tend/dt;           %total number of time steps
    p1=4;        % pN     



    %all crossbridges are initially detached:
    a=zeros(1,n0);
    
    
    for klok=1:klokmax          %loop over time steps
      t=klok*dt;                %current time
    %Note that the following ``for'' loop is commented out.
    %It explains the logic of the program but would execute slowly.
    %The code actually used is a vectorized version which appears
    %immediately below the ``for'' loop.
    % for i=1:n0                    %loop over crossbridges
    %   if(a(i))                    %if crossbridge is attached
    %     if(rand>exp(-beta*dt))    %if crossbridge detaches
    %       a(i)=0;                 %record its state as detached
    %     end
    %   else                        %crossbridge is detached
    %     if(rand>exp(-alpha*dt))         %if crossbridge attaches
    %       a(i)=1;                 %record its state as attached
    %     end
    %   end
    % end
    %The following vector statements are equivalent 
    %to the above commented-out ``for'' loop:
    
      prob=(beta*dt)*a+(alpha*dt)*(1-a);%probability of changing state
      change=-log(rand(1,n0))<prob;  %decide which crossbridges change state
      a=xor(change,a);         %change the state of those crossbridges
      
    %store results for future plotting:
      tsave(klok)=t;          %current time
      Usave(klok)=sum(a)/n0;  %fraction attached
      Psave(klok)=p1*sum(a);  %force
      asave(1:6,klok)=a(1:6);
    end
    
    
    % figure(1)
    % hold on;
    % subplot(2,1,1),stairs(tsave,Usave,'Color',C), 
    % xlabel('Time (s)'), ylabel('Fraction crossbridges attached'), hold on
    % subplot(2,1,2),stairs(tsave,Psave,'Color',C), 
    % xlabel('Time (s)'), ylabel('Force (pN)'), hold on
    % 
    % figure(2);
    % hold on;
    % pspectrum(Psave);
    % 
    % figure(3);
    % hold on;
    % subplot(length(n_totals),1,n_total);
    % histogram(Psave(tsave>.05), 2000:250:8000,'FaceColor',C) 
    % title(strcat('N=', num2str(n0)))
    % if n_total==1
    %     ylabel('Count')
    % end
    % if n_total==length(n_totals)-1
    %     xlabel('Force')
    % end
    cv = std(Psave(tsave>.05))/mean(Psave(tsave>.05));
    cov_values(n_total)= cv;
    
   
end
% subplot(5,1,2), ylabel('Crossbridge state (0=detached, 1=attached)')
% subplot(5,1,1), title('')

% 
% figure(1);
% subplot(2,1,1);
% legend(strcat('N=',string(num2cell(n_totals))))
% figure(2);
% legend(strcat('N=',string(num2cell(n_totals))))

figure();
subplot(2,1,1);
plot(n_totals, cov_values, 'b.', 'MarkerSize', 15);
xlabel('N0')
ylabel('COV (std/mean)')
subplot(2,1,2);
x = 1./sqrt(n_totals);
y = cov_values;
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(x, y, 'b.', 'MarkerSize', 15);
hold on;
plot(xFit, yFit, 'r--', 'LineWidth', 2); 
xlabel('1/sqrt(N0)')
ylabel('COV (std/mean)')