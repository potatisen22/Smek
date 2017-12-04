%Whitnoise excitation on  plexiglas measured at varius points on said glas.
%Fs, 4000Hz Calibrationfactor, 1.3 (Force/acceleration).
%chanel5 (point 13) is broken, chanel 8 has unknown point
%% A1 S12 = conj(S21)
clc;
clear all;
close all;
fs=4000;
load('fa1.mat');
x=fa1;
for i = 1:2
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
    [YF,f]=pwelch(x(:,2),WINDOW,NOVERLAP,NFFT,fs); 
    plot(f,10*log10(YF))
    hold on
    [YF,f]=pwelch(x(:,3),WINDOW,NOVERLAP,NFFT,fs);
    plot(f,10*log10(YF))
    title(['PSD estimation with block size ',num2str(NFFT), ' '])
    legend('Force','Acceleration')
    xlabel('Frequency (Hz)')
    figure()
    
    [YF,f] = cpsd(x(:,3),x(:,2),WINDOW,NOVERLAP,NFFT,fs);
    plot(f,10*log10(YF))
    hold on
    [YF,f] = cpsd(x(:,2),x(:,3),WINDOW,NOVERLAP,NFFT,fs);
    plot(f,10*log10(YF),'r') 
    title(['CrossPSD with block size',num2str(NFFT), ' '])
    legend('S_{32}','S_{23}')
end
%% A2 Calculate with FFT and welch
NFFT = 4096;
Nchunk = (2*length(x)/NFFT)-1;
beg=1;
beg2=beg;
en=NFFT;
en2=en;
WINDOW = hanning(NFFT,'periodic');
Ca = 0.5;
Cb = 3/2;
df = fs/NFFT;
NAMNARE=(3/8*df); %ca and cb thingy roughly 0.3750
f = (0:df:(fs/2)-df);
for i = 1:Nchunk %9.11 till 9.12 i boken gör vi dft per element?
    S = ((2*(abs(fft(x(beg:en,2).*WINDOW)).^2))/NAMNARE)/(NFFT^2);
    SS(:,i)=S;
    beg = beg+Nchunk;
    en=en+Nchunk;
end
Sxy=mean(SS');
figure()
semilogy(f,Sxy(1:length(Sxy)/2))
hold on
for i = 1:Nchunk %9.11 till 9.12 i boken gör vi dft per element?
    S = ((2*(abs(fft(x(beg2:en2,3).*WINDOW)).^2))/NAMNARE)/(NFFT^2);
    SS(:,i)=S;
    beg2 = beg2+Nchunk;
    en2=en2+Nchunk;
end
Sxy=mean(SS');
semilogy(f,Sxy(1:length(Sxy)/2))
title('A2, EGET PSD FÖR JAG ÄR BÄTTRE ÄN MATLAAAAAAAAAB~~~')
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
H1=(Syx./Sxx)*(1/1.3);
H2 =(Syy./Sxy)*(1/1.3);
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
legend('NFFT = 100', 'NFFT = 400')
title('LLS SINE')
window = hanning(100,'periodic');
noverlap = 100/2;
figure()
[YF,f]=pwelch(x,window,noverlap,100);
plot(f,10*log10(YF),'r');
hold on
window = hanning(400,'periodic');
noverlap = 400/2;
[YF,f]=pwelch(x,window,noverlap,400);
plot(f,10*log10(YF))
legend('NFFT = 100', 'NFFT = 400');
title('PSD SINE')
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
[YF,f]=pwelch(x,window,noverlap,4096);
plot(f,10*log10(YF),'r')
window = hanning(16384,'periodic');
noverlap = 16384/2;
hold on
[YF,f]=pwelch(x,window,noverlap,16384);
plot(f,10*log10(YF))
legend('NFFT = 4096', 'NFFT = 16384')

%% F
clc;
clear all;
close all;
fsin = 10;
t400=[1:400]'/1000;
x400=sin(2*pi*fsin*t400)+randn(size(t400))/10;
figure()
plot(t400,x400)
title('tid 400 punkter')
hold on
fsin = 10;
t450=[1:450]'/1000;
x450=sin(2*pi*fsin*t450)+randn(size(t450))/10;
plot(t450,x450)
title('tid 450 punkter')

%andra nedan
NFFT=400;
fs = 1000;
fvek=0:fs/NFFT:fs-fs/NFFT;
xfft=fft(x400,NFFT);
xfftrms = sqrt(2)*abs(xfft)/NFFT;
figure()
plot(fvek(1:end/2),xfftrms(1:end/2));
hold on
%450 below
NFFT=450;
xfft=fft(x450,NFFT);
xfftrms2 = sqrt(2)*abs(xfft)/NFFT;
plot(fvek(1:200),xfftrms2(1:200));
legend('n=400','n=450')
%% G
close all
t=[1:64000]'/1000;
x=sin(2*pi*fsin*t)+randn(size(t))/10;
%NFFT = 4096
W1=hanning(400,'periodic');
W2=hanning(450,'periodic');
x400 = fft(x(1:400),400)
x450 = fft(x(1:450),450)
x400W = fft(x(1:400).*W1,400)
x450W = fft(x(1:450).*W2,450)
fs = 1000;
fvek1=0:fs/400:fs-fs/400;
fvek2=0:fs/450:fs-fs/450;
figure()
plot(fvek1(1:length(fvek1)/2),sqrt(2)*abs(x400(1:length(fvek1)/2))/400);
hold on
plot(fvek1(1:length(fvek1)/2),sqrt(2)*abs(x400W(1:length(fvek1)/2))/400);
title ('400 blocksize')
legend('no window','window')
figure()
plot(fvek2(1:length(fvek2)/2),sqrt(2)*abs(x450(1:length(fvek2)/2))/450);
hold on
plot(fvek2(1:length(fvek2)/2),sqrt(2)*abs(x450W(1:length(fvek2)/2))/450);
title ('450 blocksize')
legend('no window','window')
%% H
close all

t=[1:64000]'/1000;
x=sin(2*pi*fsin*t)+randn(size(t))/10;
%NFFT = 4096
W1=hanning(400,'periodic');
W2=hanning(450,'periodic');
x400 = fft(x(1:400),400)
x450 = fft(x(1:450),450)
x400W = fft(x(1:400).*W1,400)
x450W = fft(x(1:450).*W2,450)
fs = 1000;
fvek1=0:fs/400:fs-fs/400;
fvek2=0:fs/450:fs-fs/450;
figure()
plot(fvek1(1:length(fvek1)/2),sqrt(2)*abs(x400(1:length(fvek1)/2))/400);
hold on
plot(fvek1(1:length(fvek1)/2),sqrt(2)*2*abs(x400W(1:length(fvek1)/2))/400);
title ('400 blocksize')
legend('no window','window')
figure()
plot(fvek2(1:length(fvek2)/2),sqrt(2)*abs(x450(1:length(fvek2)/2))/450);
hold on
plot(fvek2(1:length(fvek2)/2),sqrt(2)*2*abs(x450W(1:length(fvek2)/2))/450);
title ('450 blocksize')
legend('no window','window')
%% I1 FASFÖRSKJUTNING AV FILTER!
clc;
clear all;
close all;
f = 10;
t=[1:64000]'/1000;
x=sin(2*pi*f*t)+2* sin(3*2*pi*f*t)+ 0.5*sin(7*2*pi*f*t);
Wp = 60/(1000/2);
Ws = 70/(1000/2);
Rp = 1;
Rs = 10;
% smör
[N, Wn] = buttord(Wp, Ws, Rp, Rs);
[b,a] = butter(N,Wn);
x2 = filter(b,a,x);
fvtool(b,a);
figure()
plot(t,x2);
hold on
plot(t,x);
figure()
t2=abs(fft(x,length(x2)));
plot(t2(1:end/2))
hold on
t1 = abs(fft(x2,length(x2)));
plot(t1(1:end/2)) 
title ('Smör filter')
legend('Original signal','Damped signal')
% chernobyl
[N, Wp] = cheb1ord(Wp, Ws, Rp, Rs);
R = 0.5;
[b,a] = cheby1(N,R,Wp,'low');
x2 = filter(b,a,x);
fvtool(b,a);
figure()
plot(t,x2);
hold on
plot(t,x);
figure()
t2=abs(fft(x,length(x)));
plot(t2(1:end/2))
hold on
t1 = abs(fft(x2,length(x2)));
plot(t1(1:end/2)) 
title ('Chjernobyl filter')
legend('Original signal','Damped signal')
% Ellipserna
[N, Wp]=ellipord(Wp, Ws, Rp, Rs);
[b,a]=ellip(N,Rp,Rs,Wp,'low');
x2 = filter(b,a,x);
fvtool(b,a);
figure()
plot(t,x2);
hold on
plot(t,x);
figure()
t2=abs(fft(x,length(x)));
plot(t2(1:end/2))
hold on
t1 = abs(fft(x2,length(x2)));
plot(t1(1:end/2)) 
title ('ellip filter')
legend('Original signal','Damped signal')

%% I2 %ellips är inte så bra här egentligen, c
clc;
clear all;
close all;
load('fa1.mat')
fs = 4000;
x=fa1(:,3);
Wp = 690/(fs/2);
Ws = 710/(fs/2);
Rp = 0.1;
Rs = 20;
[N, Wp] = ellipord(Wp, Ws, Rp, Rs);
[B,A]= ellip (N,Rp,Rs,Wp,'low');
x2 = filter(B,A,x);
fvtool(B,A);
window = hanning(4096,'periodic');
noverlap = 4096/2;
figure()
[YF,f]=pwelch(x,window,noverlap,4096);
plot(f*fs/(2*pi),10*log10(YF),'r')
hold on
[YF,f]=pwelch(x2,window,noverlap,4096);
plot(f*fs/(2*pi),10*log10(YF),'b')
title ('Ellips Filter')
legend('Original signal','Damped signal')

%% J1
clc;
clear all;
close all;
f = 10;
fs = 1000
t=[1:64000]'/1000;
x=sin(2*pi*f*t)+2* sin(3*2*pi*f*t)+ 0.5*sin(7*2*pi*f*t);
Wp = [15 45]/500
Ws = [14 46]/500
Rp = 0.1;
Rs = 20;
[N, Wp] = ellipord(Wp, Ws, Rp, Rs);
[B,A]= ellip (N,Rp,Rs,Wp);
fvtool(B,A);
x2 = filter(B,A,x);
figure()
t1 = abs(fft(x2,4096));
plot(t1(1:end/2)) 
hold on
t2=abs(fft(x,4096));
plot(t2(1:end/2))
title ('ellip filter')
legend('Damped signal','Original signal')
window = hanning(4096,'periodic');
noverlap = 4096/2;
figure()
[YF,f]=pwelch(x,window,noverlap,4096);
plot(f*fs/(2*pi),10*log10(YF),'r')
hold on
[YF,f]=pwelch(x2,window,noverlap,4096);
plot(f*fs/(2*pi),10*log10(YF),'b')
title ('Ellips BANDPASS filter')
legend('Original signal','Damped signal')

%% J2
clc;
close all;
clear all;
fs = 2000;
load('fa1.mat')
x = fa1(:,3);
Wp = [500 580]/(fs/2);
Ws = [499 581]/(fs/2);
Rp = 0.1;
Rs = 20;
[N, Wp] = ellipord(Wp, Ws, Rp, Rs);
[B,A]= ellip (N,Rp,Rs,Wp);
fvtool(B,A);
x2 = filter(B,A,x);
fvtool(B,A);
window = hanning(4096,'periodic');
noverlap = 4096/2;
figure()
[YF,f]=pwelch(x,window,noverlap,4096);
plot(f*fs/(2*pi),10*log10(YF),'r')
hold on
[YF,f]=pwelch(x2,window,noverlap,4096);
plot(f*fs/(2*pi),10*log10(YF),'b')
title ('Ellips Filter')
legend('Original signal','Damped signal')