clc;
clear all;
close all;
freq=10;               % frequency of stimuli 
fs = 1024;             % Sample frequency (Hz) 
t = 0:1/fs:4-1/fs;     % 4 sec sample 
x = sqrt(2)*sin(2*pi*freq*t);
m = length(x);          % Window length 
n = pow2(nextpow2(m));  % Transform length 
y = fft(x,n);           % DFT 
figure()
plot(t,x)
figure()
df=(fs/2)/n
plot(df*(1:length(y)/2),abs(y(1:(length(y)/2)))*2)
title('onesided freq spectrum')
rmsavy=rms(x);%måste kolla i tidsdomänen
y_corr = y.*conj(y)/n^2;

plot(df*(1:length(y_corr)/2),abs(y_corr(1:(length(y_corr)/2)))*2)
title('corrected onesided freq spectrum')
%% 4.3
clc;
clear all;
close all;
k = [8, 16, 20, 32];
for i = 1:length(k)
    freq=10;               % frequency of stimuli 
    fs = k(i);             % Sample frequency (Hz) 
    t = 0:1/fs:4-1/fs;     % 4 sec sample 
    x = sqrt(2)*sin(2*pi*freq*t);
    m = length(x);          % Window length 
    n = pow2(nextpow2(m));  % Transform length 
    y = fft(x,n);           % DFT 
    figure()
    subplot(2,1,1)
    plot(t,x)
    title(['tidsdomän där ' num2str(k(i)) ' Hz'])
    df=(fs/2)/n;
    rmsavy=rms(x);%måste kolla i tidsdomänen
    y_corr = y.*conj(y)/n^2;
    subplot(2,1,2)
    plot(df*(1:length(y_corr)/2),abs(y_corr(1:(length(y_corr)/2)))*2)
    title(['corrected onesided freq spectrum ' num2str(k(i)) ' Hz'])
end