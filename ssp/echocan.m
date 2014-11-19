% ******************* Adaptive Echo Cancellation **********************
% This program reads a 16 bit mono echoed female voice sampled at
% 8000 Hz. The effect of echo on the signal was displayed and frequency
% component of the echoed signal has to be the same as the original
% signal since it is only a delayed and scaled version of the original 
% signal. This program also implements adaptive echo cancellation using 
% the LMS algortihm and produces the non-echoed signal.
%
% Filename : echocan.m  (Version 1.0)
% Programmed by Simon Hong Boon Siew
% Nanyang Technological University
% Date : 21-01-1997
% *********************************************************************

% declaration on variables, constants
	
	addech;  % call addech to obtain echo signal
      
	% Adaptive echo cancellation.

	rand('seed',123);
	L=20; b=zeros(1,L+1);

	x=echo_1+echo_2+echo_3+echo_4;	 % echo signals 
	sum_echo=x;
	d=g;        % g is wavedata with 4 echoes

	mu=0.01; alpha=0; sigma=spvari(sum_echo); px=0;

	% Adaptive LMS interference cancelling algorithm
	[y,b,px]=spnlms(x,d,b,mu,sigma,alpha,px);

	figure; 
	plot([0:N/2],y(1:N/2+1),'-g'); 
	grid;...
	title('Initial convergence of y(k)');...
	xlabel('Sample Number');...
	ylabel('Amplitude');...

	tb=[N/2+1:N];
	figure;
	plot(N/2+1:N,y(N/2+1:N),'-g'); 
	grid;...
	title('Tracking performance of y(k)');...
	xlabel('Sample Number');...
	ylabel('Amplitude');

	figure;
	diff=d-y;
	plot(t,diff);
	grid;
	sound(diff,8000);
	title('Expected output of signal after echo cancellation');...
	xlabel('Time');...
	ylabel('Amplitude');

	epsilon=max(1e-5,(diff).^2);

	figure;
	semilogy([0:N-1],epsilon,'r');...
        title('Mean Square error vs sample number');
        xlabel('Sample Number');...
        ylabel('Mean Square error')

	
	