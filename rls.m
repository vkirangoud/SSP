%/**************STATISTICAL SIGNAL PROCESSING ASSIGNMENT.*************/
%               ^^^^^^^^^^^ ^^^^^^ ^^^^^^^^^^ ^^^^^^^^^^.
% RECURSIVE LEAST SQUARES(RLS) ADAPTIVE FILTER
%1.consider the adaptive noise cancellation problem where the signal
%is given by y[n]=x[n]+v1[n],where x[n]=sin(0.05*pi*n) and v1[n]=-0.8*v1[n-1]+v[n],
%v[n] being the 0-mean,unity variance white noise sequence.consider a secondary
% source v2[n]=0.8*v2[n-1]+v[n].
%%%
% 1. simulation of the above signals.
%***********simulation of white noise sequence.************.
clc;
clear;
close all;
N = input('Enter number of iterations: '); % number of input samples.
L = input('Enter length of the filter:');      % Length of the filter.
v =  wgn(N,1,0,'real');    % 512 x 1 seqence, 0 mean,unit variance.
% power in decibels is 0 dB.
% v1(z)/v(z) = 1/(1+0.8*z^-1).
v1 = filter(1,[1,0.8],v);

% v2(z)/v(z) = 1/(1-0.8*z^-1);

v2 = filter(1,[1,-0.8],v);

% x[n] = sin(0.05*pi*n);

n = 1:N ; x(n) = sin(0.05* pi * (n-1)); % 512 samples of x[n] we have considered.
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
 V2 =flipud(buffer(v2,L,L-1)); % A buffer whose columns will act as data for each iteration. 
 % initialization of filter coefficients.
W(:,1) = zeros(L,1);
P = 0.0001 * eye(L); 
lam = 0.95; % forgetting factor.
g(1) = 0;
% Now we start iteration.
   for i = 2 : N
        g(i) = W(:,i-1)' * V2(:,i);    % output at time i.
        e(i) = y(i) - g(i);         % error signal
        K(:,i) = (P * V2(:,i)) /(lam + V2(:,i)' * P * V2(:,i));
        W(:,i) = W(:,i-1) + K(:,i) * e(i);
        P = (P - K(:,i) * V2(:,i)' * P)/lam;
    end
     figure,plot(1:N,x,1:N,e,'c--'),title('ANC-RLS');
        legend('Original','Estimated');
       z= filter(W(:,N)',1,v2');
        E = y - z;
        figure, plot(E);