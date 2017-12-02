%Whitnoise excitation on  plexiglas measured at varius points on said glas.
%Fs, 4000Hz Calibrationfactor, 1.3 (Force/acceleration).
%chanel5 (point 13) is broken, chanel 8 has unknown point
%% A1 S12 = conj(S21)
clc;
clear all;
close all;
fs=4000;
load('fa1.mat');
for i = 1:1
    x=fa1;
    if i == 1
       NFFT=4096;

    else 
        NFFT=65536;

    end
    df = fs/NFFT;
    f = (0:df:fs-df);
    NOVERLAP=NFFT/2;
    WINDOW=hanning(NFFT,'periodic');
    fs=4000;
    figure()
    pwelch(x(:,2),WINDOW,NOVERLAP,NFFT,fs); 
    hold on
    pwelch(x(:,3),WINDOW,NOVERLAP,NFFT,fs);
    figure()
    cpsd(x(:,2),x(:,3),WINDOW,NOVERLAP,NFFT,fs)
end




%% A2 Calculate with FFT and welch
close all
for i=1:2
   if i == 1
       NFFT = 4096;
   else
       
       NFFT = 65536;
   
   end
   clearvars SS
   Nchunk = (2*length(x)/NFFT)-1;
   beg=1;
   en=NFFT;
   WINDOW = hanning(NFFT,'periodic');
   Ca = 0.5;
   Cb = 3/2;
   df = fs/NFFT;
   NAMNARE=(3/8*df); %ca and cb thingy roughly 0.3750
   f = (0:df:(fs/2)-df);
   for i = 1:Nchunk %9.11 till 9.12 i boken g�r vi dft per element?
        
       S = ((2*(abs(fft(x(beg:en,2).*WINDOW)).^2))/NAMNARE)/(NFFT^2);
       SS(:,i)=S;
       beg = beg+Nchunk;
       en=en+Nchunk;
   end
   Sxy=mean(SS');
   figure()
   semilogy(f,Sxy(1:length(Sxy)/2))
end
%% B
%coherence PSD
clc;
clear all;
close all;  
load('fa1.mat');
x=fa1;
for i = 1:1
    x=fa1;
    if i == 1
       NFFT=4096;

    else 
        NFFT=65536;

    end
    NOVERLAP=NFFT/2;
    WINDOW=hanning(NFFT,'periodic');
    fs=4000;
    [Sxx,f]=pwelch(x(:,2),WINDOW,NOVERLAP,NFFT,fs);
    [Syy,f]=pwelch(x(:,3),WINDOW,NOVERLAP,NFFT,fs);
    Syx = cpsd(x(:,3),x(:,2),WINDOW,NOVERLAP,NFFT,fs);
    gammakva=(abs(Syx).^2)./(Sxx.*Syy);
    figure()
    semilogy(f,gammakva)
    hold on
end
%% C
clc;
close all; 
clear all;
load('fa1.mat');
x=fa1;
fs=4000;
NFFT = 4096;
NOVERLAP=NFFT/2;
WINDOW=hanning(NFFT,'periodic');
[Sxx,f]=pwelch(x(:,2),WINDOW,NOVERLAP,NFFT,fs);
[Syy,f]=pwelch(x(:,3),WINDOW,NOVERLAP,NFFT,fs);
Syx = cpsd(x(:,3),x(:,2),WINDOW,NOVERLAP,NFFT,fs);
Sxy =cpsd(x(:,2),x(:,3),WINDOW,NOVERLAP,NFFT,fs);
H1=Syx./Sxx;
H2 =Syy./conj(Sxy);
figure()
plot(f,abs(H1))
hold on
plot(f,abs(H2))
legend('H1','H2')
figure()
plot(f,phase(H1))
hold on
plot(f,phase(H2))

legend('phaseH1','phaseH2')
%% D
close all;
for i = 1:4
    if i == 1
       na = 4;
       nb = 3;
    elseif i == 2
       na = 8;
       nb = 6;
    elseif i == 3
       na = 16;
       nb = 14;
    else
       na = 32;
       nb = 26;
    end
    [B,A] = invfreqs(H1,2*pi*f,nb,na); 
    H=freqs(B,A,2*pi*f);
    subplot(2,2,i)
    plot(f,abs(H1),f,abs(H))
    legend('H1', 'H')
    title(['na = ', num2str(na),', nb = ',num2str(nb),])
    hold on
end
[z,p,k] = tf2zp(B,A)
figure()
plot(z,'o') 
figure()
plot(p,'x') 
%% E
%fs,fmax?
close all;
fsin = 10;
t=[1:64000]'/1000;
x=sin(2*pi*fsin*t)+randn(size(t))/10; 
for i = 1:2
    if i == 1
        NFFT=100;
    else
        NFFT=400;
    end
    xfft=fft(x,NFFT);
    sum(abs(xfft).^2);
    fs = 1000;
    fvek=0:fs/NFFT:fs-fs/NFFT;
    xfft = sqrt(2).*abs(xfft)./NFFT;
    plot(fvek(1:length(fvek)/2),abs(xfft(1:length(fvek)/2)));
    hold on

end

window = hanning(100,'periodic')
noverlap = 100/2;
figure()
pwelch(x,window,noverlap,100);
hold on
window = hanning(400,'periodic')
noverlap = 400/2;
pwelch(x,window,noverlap,400);
legend('NFFT = 100', 'NFFT = 400')
%% riktig signal
close all
load('fa1.mat');
x=fa1(:,3);
for i = 1:2
    if i == 1
        NFFT=4096;
    else
        NFFT=16384;
    end
    xfft=fft(x,NFFT);
    sum(abs(xfft).^2);
    fs = 4000;
    fvek=0:fs/NFFT:fs-fs/NFFT;
    xfft = sqrt(2).*abs(xfft)./NFFT;
    plot(fvek(1:length(fvek)/2),abs(xfft(1:length(fvek)/2)));
    hold on

end

window = hanning(4096,'periodic');
noverlap = 4096/2;
figure()
pwelch(x,window,noverlap,4096);
window = hanning(16384,'periodic');
noverlap = 16384/2;
hold on
pwelch(x,window,noverlap,16384);
legend('NFFT = 4096', 'NFFT = 16384')

%% F
