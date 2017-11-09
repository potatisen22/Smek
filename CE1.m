clc;
clear all;
close all;
%% 1.1
sound1 = wavread('sound1.wav');
sound2 = wavread('sound2.wav');
sound3 = wavread('sound3.wav');
sound4 = wavread('sound4.wav');
figure(1)
plot(sound1(:,1))
figure(2)
plot(sound2(:,1))
figure(3)
plot(sound3(:,1))
figure(4)
plot(sound4(:,1))

%% 1.2
clc;
clear all;
close all;
sound = wavread('sound1.wav');
rootmeansquare=rms(sound);
peak_value=max(sound);
ptop=peak2peak(sound);
crest=peak2rms(sound);

%% 1.3
clc;
clear all;
close all;
load 'Bearing.mat';
figure(1)
title ('Bearing 1')
histogram(Bearing1,'normalization','pdf')
figure(2)
title ('Bearing 2')
histogram(Bearing2,'normalization','pdf')
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
h=[-1, 1, 3 ,5 ,3 , 1, -1, -3,0,0,0,0];
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
