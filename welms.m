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
close all;
N = input('Enter number of iterations: '); % number of input samples.
L = input('Enter length of the filter:');      % Length of the filter.
v = wgn(N,1,0,'real');    % 512 x 1 seqence, 0 mean,unit variance.
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
 
 % Design of Length 8 - Wiener filter.
 % y[n] is the given  noisy data. x[n] is signal to be estimated from y(n).
 % v2 is input to weiner filter. v1^ is output of weiner filter. x^ = y - v1^.
 L = 8;
 y1 = v2(1:8)';
 %x1 = x(1:8)';
 % using sample autocorrelation functions.
 rv2  = (1/N) * toeplitz(y1,v2) * v2;          % v2 is 512 * 1
 rv1v2  = (1/N) * toeplitz(y1,v2) * v1;        % v1 is 512 * 1
 Rvv = toeplitz(rv2,rv2);
 h  = inv(Rvv) * rv1v2;                         % Making use of weiner Hoff equations.
 fprintf(1,'The weiner filter coefficients are: \n');
 for i = 0 : 7
     fprintf(1,'h%d = %g\n',i,h(i+1));
 end
 % plot of filtered signal , reference signal and error signal.
 xp = filter(h',1,v2');
 e = y - xp;    % estimate of x(n).
figure, subplot(2,1,1),plot(0:N-1,x',0:N-1,e,'r--');
legend('True signal','Estimated signal');
 subplot(2,1,2),plot(xp,'m');
 title('output of weiner filter');
 xlabel('sample number');

%__________________________________________________________________________________________________
%  Implementation of LMS Adaptive filter for noise cancellation.
%  --------------------------------------------------------------- 
% calculation of range of step size mu .
% ----------- -- ----- -- ---- ---- --
c = v2(1:L)';
rv2 = 1/N * toeplitz(c,v2) * v2;
evalues = eig(toeplitz(rv2));
mumax = 2/max(evalues);
fprintf(1,'\n\nenter the value of mu in the range 0 < mu < %g : \n', mumax);
mu = input('mu = ');
%---------------------------------------------------------------------

V2 =flipud(buffer(v2,L,L-1)); % A buffer whose columns will act as data for each iteration. 
[m,n] = size(V2);
% initialization of filter coefficients.
W(:,1) = zeros(L,1);
% Now we start iteration.
    for i = 1 : N
        g(i) = W(:,i)' * V2(:,i);    % output at time i.
        e(i) = y(i) - g(i);         % error signal
 
        W(:,i+1) = W(:,i) + mu * e(i) * V2(:,i) / (norm(V2(:,i),2)^2);
    end
        %
        figure,plot(1:N,x,1:N,e,'c--');
        legend('Original','Estimated');
        figure,plot(e);
        N = N+1;
        figure,subplot(2,2,1),plot(1:N,h(1),'r');
        hold on,plot(1:N,W(1,:)),title('h0 and W0');
        subplot(2,2,2),plot(1:N,h(2),'b');
        hold on,plot(1:N,W(2,:),'r'),title('h1 and W1');
        subplot(2,2,3),plot(1:N,h(3),'r');
        hold on,plot(1:N,W(3,:),'m'),title('h2 and W2');
        subplot(2,2,4),plot(1:N,h(4),'m');
        hold on,plot(1:N,W(4,:)),title('h3 and W3');
        