clc;
close all;
clear all;
[sound4,fs] = audioread('Sound4.wav');
figure();
plot(sound4(:,1));
input('hit enter to play sound')
player=audioplayer(sound4,fs);
play(player);
input('CONTINUE')
hold on
r14 = resample(sound4,1,4);
t14 = (1:length(r14)).*4;
plot(t14,r14(:,1))
hold on
r128 = resample(sound4,1,128);
t128 = (1:length(r128)).*128;
plot(t128,r128(:,1))
legend('original','1/4','1/128')
%PLOT THE ONESIDED SPECTRA HERE, DOESN'T LOOK QUITE RIGHT
figure()
df128= ((32000/2)/128)/512;%delta f
fft128=abs(fft(r128(:,1),512)); %fft
plot(df128*(1:length(fft128)/2),fft128(1:length(fft128)/2).*2);
hold on
dforg= (32000/2)/512;%delta f
fftorg = abs(fft(sound4(:,1),512)); %fft
plot(dforg*(1:length(fftorg)/2),fftorg(1:length(fftorg)/2).*2);
hold on
df14 = ((32000/2)/4)/512; %delta f
fft14 = abs(fft(r14(:,1),512)); %fft
plot(df14*(1:length(fft14)/2),fft14(1:length(fft14)/2).*2);
legend('1/128','original','1/4')
%1/128 fungerar ej
% input('1/128')
% player=audioplayer(r128,250);
% play(player);
input('1/4')
player=audioplayer(r14,8000);
play(player);