% ************************** Echo Addition ****************************
% This program reads a 16 bit mono female voice sampled at 8000Hz and 
% echos were introduced to the the signal. The effect of echo on the
% signal was displayed and the echoes are only a delayed and scaled 
% version of the original signal.	
%
% Filename : addech.m  (Version 1.0)
% Programmed by Simon Hong Boon Siew
% Nanyang Technological University
% Date : 14-12-1996
% *********************************************************************

% declaration on variables, constants
	
	load female1;
	N = length(wavedata); % wavedata size
	fs= samplingrate;
	T= N/fs;  	
	
	figure;
	t = 0 : T/(N-1) : T;
        plot(t,wavedata,'b');
	sound(wavedata,8000);
	grid;
	ylabel('Amplitude');        
	xlabel('Time domain');
     	title('16 bit data mono female speech Signal');
	
	% This program generates a signal with four echoes.
        % The original signal 'wavedata' was delayed and scaled %
	% to generate an echo signal % 
 
        echo_1 = zeros(size(t));
        echo_1(501:20800) = 0.5*wavedata(1:20300);
        echo_2 = zeros(size(t));
        echo_2(1001:20800) = 0.4*wavedata(1:19800) ;
        echo_3 = zeros(size(t));
        echo_3(5001:20800) = -0.3*wavedata(1:15800); % roll echo %
	echo_4 = zeros(size(t));
        echo_4(10001:20800) = 0.2*wavedata(1:10800);

        g = (wavedata+echo_1+echo_2+echo_3+echo_4);

	figure;
	t = 0 : T/(N-1) : T;
        plot(t,g,'b');
	sound(g,8000);
	grid;
	ylabel('Amplitude');        
	xlabel('Time domain');
     	title('Mono female speech  with 4 echo Signals');

	


