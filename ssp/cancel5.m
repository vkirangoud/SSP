% *************** Adaptive interference cancellation *****************
% This program demonstrates adaptive interference cancellation using 
% LMS algorithm. The following simulation will describe a simple
% example where the desired response d is composed of a signal s and
% a noise component nk which ideally is uncorrelated with s. The
% filter input is a noise sequence which is correlated with the noise 
% component in d. 
% 
% Filename : cancel.m  ( Version 1)
% Programmed by Simon Hong Boon Siew
% Nanyang Technological University
% Date : 22-01-1997
% *********************************************************************
 
N=501; rand('seed',123);
L=10; b=zeros(1,L+1);
mu=0.01; alpha=0; sigma=0.005; px=0;

k=[0:N-1]; arg=2*pi*k/15 + pi/3;
s=sqrt(2)*sin(2*pi*k/10) + 2*sin(2*pi*k/25);
x=0.1*sin(arg);	    % signal nk'
d=sqrt(2)*sin(2*pi*k/10) + 50*sin(2*pi*k/15) + 2*sin(2*pi*k/25) + sqrt(12)*(rand(1,N)-0.5);

figure(1);
plot([0:250],d(1:251),'y'); 
grid;...
title('Input signal dk with noise nk');...
xlabel('Time');...
ylabel('Amplitude');...

% Adaptive LMS interference cancelling algorithm
[y,b,px]=spnlms(x,d,b,mu,sigma,alpha,px);

figure(2); 
plot([0:250],y(1:251),'-w'); 
grid;...
title('Initial convergence of y(k)');...
xlabel('Time');...
ylabel('Amplitude');...

tb=[251:500];
figure(3); 
plot([251:500],y(252:501),'-w'); 
grid;...
title('Tracking performance of y(k)');...
xlabel('Time');...
ylabel('Amplitude');

figure(4);
plot(k,d-y);
grid;
title('Expected output of signal after interference cancellation');...
xlabel('Time');...
ylabel('Amplitude');

epsilon=max(1e-5,(d-y)^2);

figure(5); 
semilogy([0:N-1],epsilon,'r');...
title('Squared error vs sample number');
xlabel('Sample Number');...
ylabel('Epsilon Squared error')
