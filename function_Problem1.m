function function_Problem1
% 6880 
% Zeyu Liu
% 2/18/2020
% Adaptive filter theory 5 edition 

% Problem 1
% 1. mu = 0.01
N = 1000;
M = 10;
nums = 100;
sigma = 1;
SDV = sqrt(sigma);
mu = 0.01;
u = zeros(N,1);
e = zeros(N,1);
d = zeros(N,1);
J = zeros(N,1);
for m = 1:4
    w0 = [0.5*pi, 0.1*pi, 0.9*pi, 0.99*pi];
    for k = 1:nums %  100 independent trials of the experiment
        w = zeros(M,1);
        v = randn(N,1)*SDV;
    for n = M:N
        d(n) = 0.5*cos(w0(m)*n + pi/2) + 0.1*v(n);
        u(n) = cos(w0(m)*n);
        e(n) = d(n) - w'*u(n:-1:n-M+1);
        w = w + mu*u(n:-1:n-M+1)*e(n);
    end
    J = J + e.^2;    
    end 
    J = J/nums;

    subplot(2,2,m)
    plot(J);
    l = [0.5, 0.1, 0.9, 0.99];
    title(['w0 = ',num2str(l(m)),'\pi']);
    xlabel('Iteration');
    ylabel('Mean Squared Error');
    grid on
    grid minor
end
suptitle('Learning curve, mu = 0.01')
pause;

% 2. mu = 0.001
N = 1000;
M = 10;
nums = 100;
sigma = 1;
SDV = sqrt(sigma);
mu = 0.001;
u = zeros(N,1);
e = zeros(N,1);
d = zeros(N,1);
J = zeros(N,1);

for m = 1:4
    w0 = [0.5*pi, 0.1*pi, 0.9*pi, 0.99*pi]
    for k = 1:nums %  100 independent trials of the experiment
        w = zeros(M,1);
        v = randn(N,1)*SDV;
    for n = M:N
        d(n) = 0.5*cos(w0(m)*n + pi/2) + 0.1*v(n);
        u(n) = cos(w0(m)*n);
        e(n) = d(n) - w'*u(n:-1:n-M+1);
        w = w + mu*u(n:-1:n-M+1)*e(n);
    end
    J = J + e.^2;    
    end
    J = J/nums;
    subplot(2,2,m)
    plot(J);
    l = [0.5, 0.1, 0.9, 0.99];
    title(['w0 = ',num2str(l(m)),'\pi']);
    xlabel('Iteration');
    ylabel('Mean Squared Error');
    grid on
    grid minor
end
suptitle('Learning curve, mu = 0.001')
