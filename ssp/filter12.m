% *********************** Adaptive filter design *********************
% This program demonstrates the design of an adaptive bandpass filter
% using the LMS algorithm., The program shows that it is possible to 
% model the frequency response of a digital filter using the LMS 
% algorithm. However,the adaptive gain mu must be appropriately chosen 
% such that the adaptive process will be satisfactory. mu will depend 
% on design specifications and available computing time. 
%
% Filename : filter.m  ( Version 1)
% Programmed by Simon Hong Boon Siew
% Nanyang Technological University
% Date : 16-01-1997
% *********************************************************************

M=42; N=1001; L=100; mu=0.1; alpha=0; sigma=1+40*sqrt(2);
b=zeros(1,L+1); px=0;

amp=[zeros(1,3*(M-2)/8),0.5,ones(1,2*(M-2)/8),0.5,zeros(1,3*(M-2)/8)];
x=zeros(1,N); d=zeros(1,N);
for k=0:N-1,
   d(k+1)=sum(amp(1:42).*sin(2*pi*[0:41]*(k-L/2)/80));
   x(k+1)=sum(sin(2*pi*[0:41]*k/80));
end

[y,b,px]=spnlms(x,d,b,mu,sigma,alpha,px);
gain=abs(spgain(b,1,(0:300)*0.5/300));

figure(1); 
plot([0:300]*0.5/300,gain,'w'); 
grid;...
title('Desired and adaptive ampl. gains');
xlabel('Frequency (Hz)');...
ylabel('Amplitude gain');...
hold on;...
plot([0:41]/80,amp,'g');...
hold off

figure(2); 
stem([0:L],b,'b'); grid;...
title('Impulse response of H(z) and model');
xlabel('Sample Number');...
ylabel('Amplitude');
epsilon=max(1e-5,(d-y).^2);

figure(3); 
semilogy([0:N-1],epsilon,'r');...
title('Squared error vs sample number');
xlabel('Sample Number');...
ylabel('Epsilon squared')
