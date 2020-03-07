function function_6_19
% 6880 
% Zeyu Liu
% 2/18/2020
% Adaptive filter theory 5 edition 

% 6.19
% We need to calculate the varience of noise first
N=10000; % An enough big numbers of iteratrue will get better value
sigma=10; % guess a initial number of noise variance
SDV=sqrt(sigma);
u=zeros(N,1); 
a=-0.99; 
Tolerance=0.001; % the maximal acceptable error
MaxCycle=300; % The maximum updating cycle
error=Tolerance+1; % an initial error larger then the tolerance
mu=0.001; % like LMS learning parameter
k=1; % itialized number of cycles completed
while(abs(error)>Tolerance && k<MaxCycle)
    SDV=sqrt(sigma);   
    for n=2:N
        u(n) = -a*u(n-1)+randn(1)*SDV;
    end    
    error=(1-var(u(2:N))); % Calculate sample variance and desired variance
    sigma=sigma+mu*error; % adjust the noise variance to reduce error 
    k=k+1; % increment loop variable
end
if abs(error)>Tolerance % explaining why solution wasn't found
    fprintf('You need give a closer start value NV or more MaxCycle') 
    % systemVariance=var(u(3:N+3));
    sigma
else
    % SampleProcessVariance=var(u(3:N+3)); 
    sigma
    k
end
% sigma = 0.204
% Now we get the varience of noise sigma = 0.204
N = 300;
nums = 100;
a = -0.99;
sigma = 0.02;
SDV = sqrt(sigma);
mu = [0.01, 0.03, 0.1, 0.2, 1, 3];
u = zeros(N,1);
f = zeros(N,1);
d = zeros(N,1);
J = zeros(N,1);
for m = 1:6
    for k = 1:nums %  100 independent trials of the experiment
        w = 0;
    for n = 2:N
        u(n) = -a*u(n-1)+randn(1)*SDV;
        f(n) = u(n) - w*u(n-1);
        w = w + mu(m)*u(n-1)*f(n);
    end
    J = J + f.^2;    
    end 
    J = J/nums;
    subplot(2,3,m)
    semilogy([1:N],J);
    title(['mu = ',num2str(mu(m))]);
    xlabel('Iteration');
    ylabel('Mean Squared Error');
    grid on
    grid minor
end
suptitle('In th plot, when µ > 1 the system becomes unstable. sigma = 0.0204')
pause;