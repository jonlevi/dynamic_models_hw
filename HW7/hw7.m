
x_0 = 1:5:11;
y_0 = 1:5:11;
z_0 = 1:5:11;
t_end = 1000;
for j=1:length(x_0)
    for i=1:length(y_0)
        for q=1:length(z_0)
            y0 = [x_0(j); y_0(i); z_0(q)];
            options = odeset('RelTol',1e-8);
            
            
            [t,ysolution] = ode45(@model, [0 t_end],y0,options);
            figure(3);
            subplot(3,1,1);
            plot(t, ysolution(:,1))
            hold on;
            ylabel('W')
            subplot(3,1,2);
            
            plot(t, ysolution(:,2))
            hold on;
            ylabel('C')
            subplot(3,1,3);
            plot(t, ysolution(:,3))
            hold on;
            ylabel('F')
            xlabel('Time (min)')
            ylim([0 40])
            figure(4);
            plot3(ysolution(:,1), ysolution(:,2), ysolution(:,3));
            hold on;
            
        end
    end
end

%% 
[t2,y2] = ode45(@(t,y) 10-.5*y(1),[0 10],2,options);
figure()
plot(t2,y2);
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = model(t,y)

E=2;
d = .3;
m = .2;
b = .2;
z = .3;
n = .3;
k=.3;



dydt = [E - d*y(1)+ k*y(2)- (m*y(1)/y(2)) - b*y(1).*y(1) + n*(1-z)*y(3);
        (m*y(1)/y(2)) - k*y(2) + n*z*y(3);
        b*y(1).*y(1) - n*y(3)];

end

