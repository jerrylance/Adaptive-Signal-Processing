% ECE6880 
% Zeyu Liu
% 5/3/2020
% Final project of adaptive filter processing
% Audio Adaptive Noise Cancellation Using LMS Algorithm
function finalproject
% Read original test audio into a vector x.
[x,Fs] = audioread('original test audio.wav');
N = length(x); % N = 573300, Fs = 44100
s = x(:,1);% Extract 1 channel of the 2-channels' audio as original signal s
subplot(4,1,1);
plot(s);
sound(s,Fs);
title('Original audio s');
xlabel('Iteration');
ylabel('Amplitude');
pause;

% white noise
[n1,Fs1] = audioread('white noise.wav');
n1 = n1(:,1);
subplot(4,1,2);
plot(n1);
sound(n1,Fs1);
title('white noise n1');
xlabel('Iteration');
ylabel('Amplitude');
pause;
% pink noise
[n2,Fs2] = audioread('pink noise.wav');
n2 = n2(:,1);
subplot(4,1,3);
plot(n2);
sound(n2,Fs2);
title('pink noise n2');
xlabel('Iteration');
ylabel('Amplitude');
pause;
% car noise
[n3,Fs3] = audioread('car noise.wav');
n3 = n3(:,1);
subplot(4,1,4);
plot(n3);ylim([-1 1]);
sound(n3,Fs3);
title('car noise n3');
xlabel('Iteration');
ylabel('Amplitude');
pause;

clc;
close all;

% power spectral, LMSfilter can be used in high and low frequency signal
subplot(3,1,1)
plot((abs((fft(n1)))).^2); % white noise frequency distribution is stable
title('(c)Power Spectral of n1(white)')  
xlabel('Iteration')  
ylabel('Power Spectral Level') 
subplot(3,1,2)
plot((abs((fft(n2)))).^2); % pink noise is mainly low-frequency signals
title('(c)Power Spectral of n2(pink)')  
xlabel('Iteration')  
ylabel('Power Spectral Level') 
subplot(3,1,3)
plot((abs((fft(n3)))).^2); % car noise is mainly low-frequency signals, lower than pink noise.
title('(c)Power Spectral of n3(car)')  
xlabel('Iteration')  
ylabel('Power Spectral Level') 
pause;

clc;
close all;

% Generate the desired signal d
Fs = 44100;
d1 = s + n1; % sound with white noise
subplot(3,1,1);
plot(d1);
sound(d1,Fs);
title('sound with white noise d1');
xlabel('Iteration');
ylabel('Amplitude');
pause;
d2 = s + n2; % sound with pink noise
subplot(3,1,2);
plot(d2);
sound(d2,Fs);
title('sound with pink noise d2');
xlabel('Iteration');
ylabel('Amplitude');
pause;
d3 = s + n3; % sound with car noise
subplot(3,1,3);
plot(d3);
sound(d3,Fs);
title('sound with car noise d3');
xlabel('Iteration');
ylabel('Amplitude');
pause;
audiowrite('sound with white noise.wav',d1,Fs); % record d1.
audiowrite('sound with pink noise.wav',d2,Fs); % record d2.
audiowrite('sound with car noise.wav',d3,Fs); % record d3.

clc;
close all;

% Generate reference correlated noise n. I use white gaussian noise as reference correlated noise
rng(3); % random seed to control the reference noise
n = sqrt(0.1)*randn(N,1); % set variance = 0.1
subplot(1,1,1);
plot(n);
title('reference correlated noise n, variance = 0.1, mean = 0 white noise');
xlabel('Iteration');
ylabel('Amplitude');
pause;

clc;
close all;

% using LMS algorithm with n, d, M, mu, another function named LMSfilter.m
M = 10;% taps, filter order
mu = 0.005;% LMS learning parameter
% y is the best output sequence of LMSfilter approach correlated noise n.
[yn1,e1] = LMSfilter(n1,d1,M,mu);
[yn2,e2] = LMSfilter(n2,d2,M,mu);
[yn3,e3] = LMSfilter(n3,d3,M,mu);

% plot 
% original audio s
subplot(4,1,1);
plot(s);

sound(s,Fs);
title('Original audio s');
xlabel('Iteration');
ylabel('Amplitude');
pause;
% sound with white noise after noise reduction
subplot(4,1,2);
plot(e1);
sound(e1,Fs);
title('d1 noise reduction(white) e1');
xlabel('Iteration');
ylabel('Amplitude');
pause;
% sound with pink noise after noise reduction
subplot(4,1,3);
plot(e2);
sound(e2,Fs);
title('d2 noise reduction(pink) e2');
xlabel('Iteration');
ylabel('Amplitude');
pause;
% sound with car noise after noise reduction
subplot(4,1,4);
plot(e3);
sound(e3,Fs);
title('d3 noise reduction(car) e3');
xlabel('Iteration');
ylabel('Amplitude');
pause;
audiowrite('sound with white noise after noise redection.wav',e1,Fs); % record e1.
audiowrite('sound with pink noise after noise redection.wav',e2,Fs); % record e2.
audiowrite('sound with car noise after noise redection.wav',e3,Fs); % record e3.

clc;
close all;

% noise reduction process curve,residual noise v
v1 = s - e1;
v2 = s - e2;
v3 = s - e3;
subplot(3,1,1);
plot(v1);
title('residual noise v1');
xlabel('Iteration');
ylabel('Amplitude');
subplot(3,1,2);
plot(v2);
title('residual noise v2');
xlabel('Iteration');
ylabel('Amplitude');
subplot(3,1,3);
plot(v3);
title('residual noise v3');
xlabel('Iteration');
ylabel('Amplitude');
% sound(v1,Fs); % we don't care about the sound of residual noise.
pause;

clc;
close all;

% Sampling comparison noise reduction
t = (400000:400500);
subplot(3,1,1);
plot(t,e1(400000:400500),'k',t,v1(400000:400500),t,s(400000:400500),'g','linewidth',1);grid;
title('Sampling comparison of noise(white) reduction');
legend('e1,sound after noise reduction','v1,residual noise','s,original audio','Location','southeast')
xlabel('Iteration');
ylabel('Amplitude');
subplot(3,1,2);
plot(t,e2(400000:400500),'k',t,v2(400000:400500),t,s(400000:400500),'g','linewidth',1);grid;
title('Sampling comparison of noise(pink) reduction');
legend('e2,sound after noise reduction','v2,residual noise','s,original audio','Location','southeast')
xlabel('Iteration');
ylabel('Amplitude');
subplot(3,1,3);
plot(t,e3(400000:400500),'k',t,v3(400000:400500),t,s(400000:400500),'g','linewidth',1);grid;
title('Sampling comparison of noise(white) reduction');
legend('e3,sound after noise reduction','v3,residual noise','s,original audio','Location','southeast')
xlabel('Iteration');
ylabel('Amplitude');
pause;

clc;
close all;

% Use the best estimate weight to get the best output sequence yn
% yn is approximately equal to reference correlated noise n
subplot(3,1,1);
plot(yn1);
title('yn1,the best output sequence approximately equal to reference noise n(white)');
xlabel('Iteration');
ylabel('Amplitude');
subplot(3,1,2);
plot(yn2);
title('yn2,the best output sequence approximately equal to reference noise n(pink)');
xlabel('Iteration');
ylabel('Amplitude');
subplot(3,1,3);
plot(yn3);
title('yn3,the best output sequence approximately equal to reference noise n(car)');
xlabel('Iteration');
ylabel('Amplitude');
pause;

clc;
close all;

% To prove that the stronger the correlation, the better the noise reduction effect.
% I generate AR model noise as signal noise and MA model noise as reference correlated noise n, 
% we can adjust AR and MA parameters to improve the result.

% Generating noise of AR model
ar = [1, 1/2];   % AR parameter
n_ar=filter(1,ar,n);
subplot(2,1,1);
plot(n_ar);
title('AR model noise, n_a_r');
xlabel('Iteration');
ylabel('Amplitude');
% Generating noise for MA model, it's the correlated noise of AR model noise
ma = [1,-0.8,0.4,-0.2];   % MA parameter
n_ma = filter(ma,1,n);
subplot(2,1,2);
plot(n_ma);
title('MA model noise, reference correlated noise, n_m_a');
xlabel('Iteration');
ylabel('Amplitude');
pause;

clc;
close all;

Fs = 44100;
da = s + n_ar; % sound with AR model noise
subplot(3,1,1);
plot(da);
sound(da,Fs);
audiowrite('sound with AR model noise.wav',da,Fs); % record da.
title('sound with AR model noise, da');
xlabel('Iteration');
ylabel('Amplitude');
pause;
subplot(3,1,2);
[~,ea] = LMSfilter(n_ma,da,M,mu);
plot(ea);
sound(ea,Fs); % sound after noise reduction
audiowrite('sound with AR model noise after noise reduction.wav',ea,Fs); % record ea.
title('da noise reduction(AR), ea');
xlabel('Iteration');
ylabel('Amplitude');
pause;
va = s - ea; % residual noise va
subplot(3,1,3);
plot(va);
title('residual noise va');
xlabel('Iteration');
ylabel('Amplitude');
% SNR
r1 = snr(da,n_ar);
r2 = snr(ea,va);
disp(r1); % r1 = 0.3929
disp(r2); % r2 = 17.9667
pause;

clc;
close all;

% When change mu or M, the result also changed. Using AR-MA model noise.
% M is the same, change mu
M = 10;
mu = [0.001,0.005,0.1];
for m = 1:3
    subplot(3,2,2*m-1);
    [~,emu] = LMSfilter(n_ma,da,M,mu(m));
    plot(emu);ylim([-1.2 1.2]);
    title(['after noise reduction,emu, \mu = ',num2str(mu(m))]);
    xlabel('Iteration');
    ylabel('Amplitude');
    vmu = s - emu; % residual noise vmu
    subplot(3,2,2*m);
    plot(vmu);
    title(['residual noise,vmu, \mu = ',num2str(mu(m))]);
    xlabel('Iteration');
    ylabel('Amplitude');
    sound(emu,Fs); % sound after noise reduction
    pause;
end

clc;
close all;

% mu is the same, change M
M = [5,10,30];
mu = 0.005;
for m = 1:3
    subplot(3,2,2*m-1);
    [~,emu1] = LMSfilter(n_ma,da,M(m),mu);
    plot(emu1);ylim([-1.2 1.2]);
    title(['after noise reduction,emu, M = ',num2str(M(m))]);
    xlabel('Iteration');
    ylabel('Amplitude');
    vmu1 = s - emu1; % residual noise vmu
    subplot(3,2,2*m);
    plot(vmu1);
    title(['residual noise,vmu, M = ',num2str(M(m))]);
    xlabel('Iteration');
    ylabel('Amplitude');
    sound(emu1,Fs); % sound after noise reduction
    pause;
end

clc;
close all;
% learning curve J, find the MSE of different mu and M.
% However, It's too slow to calculate all MSE, so just calculate some point MSE.
% M is the same, change mu
M = 10;
mu = [0.001,0.005,0.1];
h = [200,15000,50000]; % change h to calculate the MSE in h iteration.
for mm = 1:3
subplot(3,2,2*mm-1);    
    for m = 1:3
        [~,emu] = LMSfilter(n_ma,da,M,mu(m));
        Ja = sum((s - emu).^2);
        J = ((s - emu).^2 + Ja)./h(mm);
        plot(J(1:h(mm)));
        mean(J)
        hold on
    end
title('MSE') 
xlabel(['When Iteration = ',num2str(h(mm))]);
ylabel('Mean Square Error');
legend('mu = 0.001','mu = 0.005','mu = 0.1')
hold off
end
% mu is the same, change M
M = [5,10,30];
mu = 0.005;
for mm = 1:3
subplot(3,2,2*mm);
    for m = 1:3
        [~,emu] = LMSfilter(n_ma,da,M(m),mu);
        Ja = sum((s - emu).^2);
        J = ((s - emu).^2 + Ja)./h(mm);
        plot(J(1:h(mm)));
        mean(J)
        hold on
    end
title('MSE') 
xlabel(['When Iteration = ',num2str(h(mm))]);
ylabel('Mean Square Error')
legend('M = 5','M = 10','M = 30')
hold off
end
pause;
% figure1: 2.7945, 0.5892, 1.0023
% figure2: 1.2461, 0.5892, 0.6349
% figure3: 0.0373, 0.0079, 0.0134
% figure4: 0.0166, 0.0079, 0.0085
% figure5: 0.0112, 0.0024, 0.0040
% figure6: 0.0050, 0.0024, 0.0025
