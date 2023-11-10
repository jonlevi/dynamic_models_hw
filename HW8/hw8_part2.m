%% problem 2

% uniform distribution
N = 100000;
r = 2.*rand(N,4) - 1;
N_saddles = 0;
N_stable_nodes = 0;
N_unstable_nodes = 0;
N_stable_spirals = 0;
N_unstable_spirals = 0;
N_eig_zero = 0;
error=0;

for i=1:N
    A = [r(i,1) r(i,2); r(i,3) r(i,4)];
    lambdas = eigs(A);
    l1 = lambdas(1);
    l2 = lambdas(2);
    if (l1==0) || (l2==0)
        N_eig_zero = N_eig_zero + 1;
        continue
    end

    isReal = (isreal(l1)) && (isreal(l2));
    isStable = (real(l1) < 0) && (real(l2)<0);
    isUnstable = (real(l1) > 0) && (real(l2)>0);
    isSaddle = isReal && (real(l1)*real(l2) < 0);
    if isSaddle
        N_saddles = N_saddles + 1;
        continue
    end
    
    if isReal
        if isStable
            N_stable_nodes = N_stable_nodes+1;
        elseif isUnstable
            N_unstable_nodes = N_unstable_nodes+1;
        else
            error = error+1;
        end
        
        continue

    else
        if isStable
            N_stable_spirals = N_stable_spirals+1;
        elseif isUnstable
            N_unstable_spirals = N_unstable_spirals+1;
        else
            error = error+1;
        end
        
        continue

    end

end

% figure();
x = ["Saddle" "Stable Node" "Unstable Node" "Stable Spiral" "Unstable Spiral"];
y1 = [N_saddles N_stable_nodes N_unstable_nodes N_stable_spirals N_unstable_spirals];
% bar(x,y./N)
% ylabel('Frequency','FontSize',16)
% xlabel('Type of Fixed Point', 'FontSize',16)
% title('Stability of Uniform 2D Random System','FontSize',16)

%%
% normal distribution
N = 100000;
r = randn(N,4);
N_saddles = 0;
N_stable_nodes = 0;
N_unstable_nodes = 0;
N_stable_spirals = 0;
N_unstable_spirals = 0;
N_eig_zero = 0;
error=0;

for i=1:N
    A = [r(i,1) r(i,2); r(i,3) r(i,4)];
    lambdas = eigs(A);
    l1 = lambdas(1);
    l2 = lambdas(2);
    if (l1==0) || (l2==0)
        N_eig_zero = N_eig_zero + 1;
        continue
    end

    isReal = (isreal(l1)) && (isreal(l2));
    isStable = (real(l1) < 0) && (real(l2)<0);
    isUnstable = (real(l1) > 0) && (real(l2)>0);
    isSaddle = isReal && (real(l1)*real(l2) < 0);
    if isSaddle
        N_saddles = N_saddles + 1;
        continue
    end
    
    if isReal
        if isStable
            N_stable_nodes = N_stable_nodes+1;
        elseif isUnstable
            N_unstable_nodes = N_unstable_nodes+1;
        else
            error = error+1;
        end
        
        continue

    else
        if isStable
            N_stable_spirals = N_stable_spirals+1;
        elseif isUnstable
            N_unstable_spirals = N_unstable_spirals+1;
        else
            error = error+1;
        end
        
        continue

    end

end

figure();
x = ["Saddle" "Stable Node" "Unstable Node" "Stable Spiral" "Unstable Spiral"];
y2 = [N_saddles N_stable_nodes N_unstable_nodes N_stable_spirals N_unstable_spirals];
bar(x,[y1./N; y2./N]);
ylabel('Frequency','FontSize',16)
xlabel('Type of Fixed Point', 'FontSize',16)
title('Stability of Random 2D System','FontSize',16)

a = legend('$$U[-1, 1]$$' ,'$$\eta(0,1)$$');
set(a,'Interpreter','latex','fontsize',14)


