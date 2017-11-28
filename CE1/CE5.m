clc;
clear all;
close all;
%sampling intervall dt = 0.2 s pixle size 1mm
%5.1 autocorr
im1=rgb2gray(imread('Correlation_2D_1.jpg'));
im2=rgb2gray(imread('Correlation_2D_2.jpg'));
%% 5.1
figure()
mesh(xcorr2(im1));
title('autocorr pic 1')
figure()
mesh(xcorr2(im2));
title('autocorr pic 2')
%% 5.2
xcor=xcorr2(im1,im2);
mesh(xcor);
t = 0;
for i = 1:length(xcor(:,1))
    for k =1:length(xcor(1,:))
        if xcor(i,k) > t
            t=xcor(i,k);
            ii=i;
            kk=k;
        end
    end
end
lengtravx=(abs(ii)-270);
lengtravy=(abs(kk)-346);
sx = (lengtravx/0.2)/1000;
sy = (lengtravy/0.2)/1000;
