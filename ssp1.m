%/**************STATISTICAL SIGNAL PROCESSING ASSIGNMENT.*************/
%               ^^^^^^^^^^^ ^^^^^^ ^^^^^^^^^^ ^^^^^^^^^^.
%1.consider the adaptive noise cancellation problem where the signal
%is given by y[n]=x[n]+v1[n],where x[n]=sin(0.05*pi*n) and v1[n]=-0.8*v1[n-1]+v[n],
%v[n] being the 0-mean,unity variance white noise sequence.consider a secondary
% source v2[n]=0.8*v2[n-1]+v[n].
%%%
% 1. simulation of the above signals.
%***********simulation of white noise sequence.************.
clc;
clear;
v = wgn(512,1,0,'real');    % 512 x 1 seqence, 0 mean,unit variance.
% power in decibels is 0 dB.
% v1(z)/v(z) = 1/(1+0.8*z^-1).
v1 = filter(1,[1,0.8],v);

% v2(z)/v(z) = 1/(1-0.8*z^-1);

v2 = filter(1,[1,-0.8],v);

% x[n] = sin(0.05*pi*n);

n = 1:512 ; x(n) = sin(0.05* pi * (n-1)); % 512 samples of x[n] we have considered.
 y(n) = x(n)' + v1(n);        % x[n]' will have same dimensions as that of v1[n].
 
 subplot(5,1,1),plot(v);text(514.9874,3.4659,'V[n]');
 subplot(5,1,2);
 grid on;
 plot(v1,'r');text(530.1008,1.2931,'V1[n]');
 subplot(5,1,3),plot(x,'m');text(511.9647,0.7356,'X[n]=sin(0.05pi*n)');
 subplot(5,1,4),plot(y,'b');text(515.7431,6.6379,'Y[n]');
 subplot(5,1,5),plot(v2,'g');text(516.4987,2.0690,'V2[n]');
 %legend('xaxis-sample number','yaxis-amplitude');
 xlabel('samples');
 %*****************************************************************************
 %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
 % Design of Length 8 - Wiener filter.
 % y[n] is the given  noisy data. x[n] is reference signal.
 y1 = y(1:8)';
 x1 = x(1:8)';
 % using sample autocorrelation functions.
 ryy  = (1/512) * toeplitz(y1,y) * y';          % y' becomes 512 * 1
 rxy  = (1/512) * toeplitz(y1,y) * x';        % x' becomes 512 * 1
 Ryy = toeplitz(ryy,ryy);
 h  = inv(Ryy) * rxy;
 fprintf(1,'The weiner filter coefficients are: \n');
 disp(h');
 % plot of filtered signal , reference signal and error signal.
 xp = filter(h',1,y');
 e = x' - xp;
figure, subplot(2,1,1),plot(0:511,x',0:511,xp,'r--');
legend('reference signal','estimated signal');
 subplot(2,1,2),plot(e,'m');
 title('error signal');
 xlabel('sample number');
 
 