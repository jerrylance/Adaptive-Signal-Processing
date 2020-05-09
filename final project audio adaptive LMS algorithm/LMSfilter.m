% ECE6880 
% Zeyu Liu
% 5/3/2020
% Final project of adaptive filter processing
% Audio Adaptive Noise Cancellation Using LMS Algorithm
function [yn,e]=LMSfilter(xn,d,M,mu)
N = length(xn);
e = zeros(N,1); % error between prediction and true results
W = zeros(M,N); % M:taps, W:Filter weight matrix, Initial value is 0

for n = M:N     % n-th iteration
    x = xn(n:-1:n-M+1); % M taps inputs
    y = W(:,n-1)'*x; % output of LMS filter
    e(n) = d(n)-y; % n-th error of iteration
    W(:,n)=W(:,n-1)+mu*x*e(n); % update weight
end

% The best output sequence of the optimal filter
yn = ones(size(xn)); 
for k = M:N
    xb = xn(k:-1:k-M+1);
    yn(k) = W(:,end)'* xb; % Use the best estimate weight get the output
end
