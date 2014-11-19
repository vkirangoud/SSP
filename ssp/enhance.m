% *********************** Adaptive line enhancement *******************
% This program demostrates the concept of an adaptive line enhancemnt
% structure. It adapts an input signal corrupted by white noise and for
% this white noise, a delay of 1 unit is sufficient to decorrelate the 
% noise in dk and xk. In addition, we will demonstrate the effect of
% changing the convergence parameter u on the adaptive line enhancement 
% of the noisy input sequence dk.
%
% Filename : enhance.m  ( Version 1)
% Programmed by Simon Hong Boon Siew
% Nanyang Technological University
% Date : 12-01-1997
% *********************************************************************

mu=0.1; sigma=2; alpha=0; L=20; N=501;
b=zeros(1,L+1); px=0;
% input noisy sequence with white noise
d=sqrt(2)*sin(2*pi*[0:N-1]/20) + sqrt(12)*(rand(1,N)-0.5);

% x is input delay version of d
x=[0,d(1:N-1)];

[y,b,px]=spnlms(x,d,b,mu,sigma,alpha,px);

figure(1);
plot([0:500],d,'y');
grid;...
title('Noisy input sequence, d(k)');...
xlabel('Time');...
ylabel('Amplitude');...

figure(2);
plot([0:500],y,'g'); 
grid;...
title('Output sequence y(k); mu=0.1');...
xlabel('Time');...
ylabel('Amplitude');

b=zeros(1,L+1); 
px=0; mu=0.01;
[y,b,px]=spnlms(x,d,b,mu,sigma,alpha,px);

figure(3);
plot([0:500],y,'r'); 
grid;...
title('Output sequence y(k); mu=0.01');...
xlabel('Time');...
ylabel('Amplitude');
