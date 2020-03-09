function function_7_10
% 6880 
% Zeyu Liu
% 2/25/2020
% Adaptive filter theory 5 edition 

% 7.10
% (a)(b)(c)
N = 200;
nums = 100; % Montecarlo
a1 = 0.1; % AR parameter
a2 = -0.8;
NV = 0.2775; % Noise variance from the 6.17(a),perivous homework5. So ignored code.
SDV = sqrt(NV);
mu = 0.2;
delta = [0.5,0.25,0.75];

u = zeros(N+3,1);
f = zeros(N+3,1); 
e1 = zeros(N+3,1);
e2 = zeros(N+3,1);
J = zeros(N+3,1);  

for m = 1:3
    s = rng(3);% give a random seed
    for k = 1:nums
        W = zeros(2,N+3);
        for n = 3:N+3 
            u(n) = -a1*u(n-1)-a2*u(n-2)+randn(1)*SDV;
            f(n) = u(n)-W(1,n-1)*u(n-1)-W(2,n-1)*u(n-2);
            e1(n) = a1-W(1,n-1);
            e2(n) = a2-W(2,n-1);
            W(:,n)=W(:,n-1)+(mu/(delta(m)+u(n-1)^2+u(n-2)^2))*[u(n-1);u(n-2)]*f(n);
        end
        J = J+f.^2;
    end
    J = J/nums;
    
    x = 1:N+3;
    subplot(3,3,1+(m-1)*3)
    plot(x,J,'linewidth',1);xlim([0 N+3]);ylim([0 0.6])
    title(['NLMS Learning curve \delta = ',num2str(delta(m))])
    xlabel('Iterations') 
    ylabel('Mean Squared Error')
    grid on
    grid minor
    subplot(3,3,2+(m-1)*3)
    plot(x,e1,'linewidth',1,'color','#BF360C');xlim([0 N+3]);ylim([-0.2 0.6])
    title(['\epsilon_1 Tap-weight estiamte Error \delta = ',num2str(delta(m))]) 
    xlabel('Iterations') 
    ylabel('error between weight and a1 \epsilon_1')
    grid on
    grid minor
    subplot(3,3,3+(m-1)*3)
    plot(x,e2,'linewidth',1,'color','#3E2723');xlim([0 N+3]);ylim([-2 0])
    title(['\epsilon_2 Tap-weight estiamte Error \delta = ',num2str(delta(m))]) 
    xlabel('Iterations') 
    ylabel('error between weight and a2 \epsilon_2')
    grid on
    grid minor
end
% when ? increases, the error curve becomes more stable, the robustness grows.
% But the learning curve become more unstable, the fluctuation grows.