function function_6_17
% 6880 
% Zeyu Liu
% 2/18/2020
% Adaptive filter theory 5 edition 

% 6.17
% (a)
N=10000; % An enough big numbers of iteratrue will get better value
NV=10; % guess a initial number of noise variance
SDV=sqrt(NV);
u=zeros(N+3,1); 
a1=0.1; % AR parameter
a2=-0.8;
Tolerance=0.001; % the maximal acceptable error
MaxCycle=300; % The maximum updating cycle
error=Tolerance+1; % an initial error larger then the tolerance
mu=0.01; % like LMS learning parameter
k=1; % itialized number of cycles completed
while(abs(error)>Tolerance && k<MaxCycle)
    SDV=sqrt(NV);   
    for n=3:(N+3)
        u(n)=-a1*u(n-1)-a2*u(n-2)+randn(1)*SDV;
    end    
    error=(1-var(u(3:N+3))); % Calculate sample variance and desired variance
    NV=NV+mu*error; % adjust the noise variance to reduce error 
    k=k+1; % increment loop variable
end
if abs(error)>Tolerance % explaining why solution wasn't found
    fprintf('You need give a closer start value NV or more MaxCycle') 
    % systemVariance=var(u(3:N+3));
    NV
else
    % SampleProcessVariance=var(u(3:N+3)); 
    NV
    k
end
% NV = 0.2775

% (b)
N = 200;
nums = 100;
a1 = 0.1; % AR parameter
a2 = -0.8;
NV = 0.2775;
SDV = sqrt(NV);
mu = 0.05;

u = zeros(N+3,1);
f = zeros(N+3,1); 
A1 = zeros(N+3,1);
A2 = zeros(N+2,1);
J = zeros(N+3,1); 

for k = 1:nums
    W = zeros(2,N+3);
    for n = 3:N+3
        u(n) = a1*u(n-1)+a2*u(n-2)+randn(1)*SDV;
        f(n) = u(n)-W(1,n-1)*u(n-1)-W(2,n-1)*u(n-2);
        A1(n) =W(1,n-1);
        A2(n) =W(2,n-1);
        W(:,n)=W(:,n-1)+mu*f(n)*[u(n-1);u(n-2)];
    end
    J = J+f.^2;
end
st1 = sprintf('The estimated a1 is %.4f\n', sum(A1(100:N+3))/104); 
st2 = sprintf('The estimated a2 is %.4f\n', sum(A2(100:N+3))/104);
J = J/nums; 
x = 1:N+3;
subplot(2,2,[1,2])
plot(x,J,'linewidth',1);grid on;grid minor;xlim([0 N+3])
title('Learning curve, noisevariance = 0.2775');xlabel('Iterations');ylabel('Mean Squared Error')
subplot(2,2,3)
plot(x,A1,'linewidth',1);grid on;grid minor;xlim([0 N+3])
title('Process of estimate a1');xlabel(st1);ylabel('estiamte a1')
subplot(2,2,4)
plot(x,A2,'linewidth',1);grid on;grid minor;xlim([0 N+3])
title('Process of estimate a2');xlabel(st2);ylabel('estiamte a2')
pause;

% (c)(d)(e)
N = 200;
nums = 100;
a1 = 0.1; % AR parameter
a2 = -0.8;
NV = 0.2775;
SDV = sqrt(NV);
mu = 0.05;

u = zeros(N+3,1);
f = zeros(N+3,1); 
e1 = zeros(N+3,1);
e2 = zeros(N+3,1);
J = zeros(N+3,1);  % experiment error
JJ = zeros(N+3,1); % theoretical error

for k = 1:nums
    W = zeros(2,N+3);
    for n = 3:N+3
        u(n) = a1*u(n-1)+a2*u(n-2)+randn(1)*SDV;
        f(n) = u(n)-W(1,n-1)*u(n-1)-W(2,n-1)*u(n-2);
        e1(n) = a1-W(1,n-1);
        e2(n) = a2-W(2,n-1);
        W(:,n)=W(:,n-1)+mu*f(n)*[u(n-1);u(n-2)];
    end
    J = J+f.^2;
end
J = J/nums; % experiment error
for n = 3:N+3
    JJ(n) = (1-NV*(1+mu/2))*(1-mu)^(2*n)+NV*(1+mu/2) % theoretical error
end
x = 1:N+3;
grid on
grid minor
subplot(2,2,[1,2])
plot(x,(abs((fft(f)))).^2,'linewidth',1);xlim([0 N+3])
title('(c)Power Spectral of f') 
xlabel('Iterations') 
ylabel('Power Spectral Level')
grid on
grid minor
subplot(2,2,3)
plot(x,(abs((fft(e1)))).^2,'linewidth',1);xlim([0 N+3])
title('Power Spectral of \epsilon_1') 
xlabel('Frequency') 
ylabel('Power Spectral Level')
grid on
grid minor
subplot(2,2,4)
plot(x,(abs((fft(e2)))).^2,'linewidth',1);xlim([0 N+3])
title('Power Spectral of \epsilon_2') 
xlabel('Frequency') 
ylabel('Power Spectral Level')
grid on
grid minor
pause;

%(e)
subplot(1,1,1)
plot(x,J,x,JJ,'linewidth',1);xlim([0 N+3])
legend('Experimental Error','Theoretical Error')
title('Compare With Experimental And Theoretical Result') 
xlabel('Iteration') 
ylabel('Mean Square Error')
grid on
grid minor



