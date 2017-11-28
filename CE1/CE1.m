%1.1
clc;
clear all;
close all;
[y1,fs1] = audioread('sound1.wav');
[y2,fs2] = audioread('sound2.wav');
[y3,fs3] = audioread('sound3.wav');
[y4,fs4] = audioread('sound4.wav');
%använder int16 då negativa tal alltså 16bits per sample, 2 kanaler!!!!
bits1 = length(y1)*16*2;
bits2 = length(y2)*16*2;
bits3 = length(y3)*16*2;
bits4 = length(y4)*16*2;
figure(1)
plot(y1)
figure(2)
plot(y2)
figure(3)
plot(y3)
figure(4)
plot(y4)

%% 1.2
clc;
clear all;
close all;
sound = audioread('sound1.wav');
rootmeansquare=rms(sound);
peak_value=max(abs(sound));
ptop=peak2peak(sound);
crest=peak2rms(sound);

%% 1.3
clc;
clear all;
close all;
load 'Bearing.mat';
figure(1)
histogram(Bearing1,'normalization','pdf')
title ('Bearing 1 probability density function')


figure(2)
histogram(Bearing2,'normalization','pdf')
title ('Bearing 2 probability density function')


[n1,x1]=histcounts(Bearing1,'normalization','pdf');
dx1=abs(max(x1)-min(x1))/(length(x1)-1);
sumpdf1=sum(n1*dx1);
[n2,x2]=histcounts(Bearing2,'normalization','pdf');
dx2=abs(max(x2)-min(x2))/(length(x2)-1);
sumpdf2=sum(n2*dx2);
if sumpdf2 == sumpdf1
    print('bearing 1 och 2 har båda Fördelningsfunktion med summan 1')
else
    print('ACHTUNG ER IST EIN PROBLEM IN DER KOD DIE SUM IM DERFÖRDELNINGSFUNKTION IST NICHT EIN')
end
M1 = max(n1)
M2 = max(n2)
%skewness och kurtosis below
sigmab1 = std(Bearing1);
sigmab2 = std(Bearing2);
N1=length(Bearing1);
N2=length(Bearing2);
M3b1 =1/N1*sum(Bearing1.^3);
M4b1 =1/N1*sum(Bearing1.^4);
M4b2 =1/N2*sum(Bearing2.^4);
M3b2 =1/N2*sum(Bearing2.^3);
skewb1=M3b1/sigmab1^3;
qurtb1=M4b1/sigmab1^4;
skewb2=M3b2/sigmab2^3;
qurtb2=M4b2/sigmab2^4;
b = [skewb1,qurtb1 ; skewb2 qurtb2]

%% 1.4
clear all;
close all;
clc;
x1org =[0,1,0,0,0,0]
horg = [-1,1,3,4,5,1,-1,-3]
x1=[0,1,0,0,0,0,0,0,0,0,0,0,0];
%h=[-1, 1, 3 ,5 ,3 , 1, -1, -3,0,0,0,0];
horg = [-4,6,2,0,0,0,3,0.1,0,0,0,0]
h = horg
for n = 1:length(h)
   for i=1:n
      y1(i)= x1(i)*h(n-i+1); 
   end
   
   y1con(n)=sum(y1);
end
y1 = conv(x1org,horg);
y1con;
x2 = [0,1,0,0,1,0];
x3 =2*sin((0:0.1:2*pi));
y2 = conv(x2,horg);
y3 = conv(x3,horg);
figure(1)
stem(y1)
title('y with impulse x1')
figure(2)
stem(y2)
title('y with impulse x2')
figure(3)
stem(y3)
title('y with impulse x3')
%% generera h-vektor
dt=0.001;
n=5;
for i = 1:201
    t = dt*(i-1)
    htruck(i) =sin(100*t)*exp(-n*t);
end
xtruck = zeros(1,101)
xtruck(101)=1
figure(4)
stem(conv(xtruck,htruck))
title('(output) truck n=5')
figure(5)
stem(htruck)
title('system response fucntion n=5')
figure(6)
stem(xtruck)
title('impulse function n = 5')

%% n=50 nedan
n=50;
for i = 1:201
    t = dt*(i-1);
    htruck(i) =sin(100*t)*exp(-n*t);
end
xtruck = zeros(1,101)
xtruck(101)=1
figure(7)
stem(conv(xtruck,htruck))
title('truck n=50')
figure(8)
stem(htruck)
title('system response fucntion n=50')
figure(9)
stem(xtruck)
title('impulse function n = 50')