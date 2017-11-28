%% 2.1 a och b
clc;
close all;
clear all;
time=0.001*[1:1000]';
sin_4_harmonic=1*sin(time*4*pi)+0.25*sin(time*4*4*pi);
[cN1,gN1]=fit(time, sin_4_harmonic,'fourier1');
[cN2,gN2]=fit(time, sin_4_harmonic,'fourier2');
[cN4,gN4]=fit(time, sin_4_harmonic,'fourier4');
[cN8,gN8]=fit(time, sin_4_harmonic,'fourier8');
figure(1)
plot(cN1,'r')
hold on
plot(cN2,'y')
hold on
plot(cN4,'g')
hold on
plot(cN8,'b')
hold on
plot(time,sin_4_harmonic)
legend('N=1','N=2','N=4','N=8','original sinusfunktion')

[y,f] = audioread('Sound1.wav');
t=1:length(y);
[cN1,gN1]=fit(t',y(:,1),'fourier1');
[cN2,gN2]=fit(t',y(:,1),'fourier2');
[cN4,gN4]=fit(t',y(:,1),'fourier4');
[cN8,gN8]=fit(t',y(:,1),'fourier8');
figure(2)
plot(t,y(:,1),'y')
hold on
plot(cN1,'r')
hold on
plot(cN2,'k')
hold on
plot(cN4,'g')
hold on
plot(cN8,'b')
legend('original funktion','N=1','N=2','N=4','N=8')
%% 2.2
clc;
close all;
clear all;
[y,f] = audioread('Sound1.wav');
df = 32000*(length(y)/f)/(512*2)/(length(y)/f);
df500=df;
L = 500;
NFFT = 2^nextpow2(L);
k = fft(y(:,1),NFFT);
sk=fftshift(k);

t = -256:255;
figure(1)
plot(t.*df,abs(sk));
title('L=500')
%L = 1000
L = 1000;
df = 32000*(length(y)/f)/(1024*2)/(length(y)/f);
df1000=df;
NFFT = 2^nextpow2(L);
k1000 = fft(y(:,1),NFFT);
sk1000=fftshift(k1000);
figure(2)
t = -512:511;
plot(t.*df,abs(sk1000));
title('L=1000')
%ifft
ift500 = ifft(k);
ift1000 = ifft(k1000);
figure()
plot(y(:,1))
hold on
plot(ift1000)
hold on
plot(ift500) %SKITEN SKÄRS JU BORT!!!!!!!
%% 2.3
clc;
close all;
clear all;
[y,f] = audioread('Sound1.wav');
figure(1)
spectrogram(y(:,1),64,0,512,f);
title('64')

figure(2)
spectrogram(y(:,1),8192,0,512,f);

title('8192')