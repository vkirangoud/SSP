% % /**************STATISTICAL SIGNAL PROCESSING ASSIGNMENT.*************/
%               ^^^^^^^^^^^ ^^^^^^ ^^^^^^^^^^ ^^^^^^^^^^.
%2. Channel Equalisation.
% h[n] = 0.1 * (0.5)^n, n =0,1,...8.
%x training signal.
clear;
clc;
close all;
N = 200;
% Generate bit stream

x = rand(N,1);
z0 = find(x < 0.5);
z1 = find(x >= 0.5);
x(z0) = -1*ones(size(z0));
x(z1) = +1*ones(size(z1));
    
n=1:9; h = 0.1 * ((0.5) .^ (n-1));
y = filter(h,1,x);
% yn = y + sqrt(0.005) * randn(1,N);
yn = y + normrnd(0,0.005,N,1);
L = 20;
% mu = 0.15;
V2 = flipud(buffer(yn,L,L-1)); % A buffer whose columns will act as data for each iteration
W(:,1) = zeros(L,1);
P = 0.0001 * eye(L); 
V2 = [zeros(L,1) V2];
lam = 0.95; % forgetting factor.
% Now we start iteration.
    g(1) = 0;
   for i = 1 : N-1
        g(i+1) = W(:,i)' * V2(:,i+1);    % output at time i.
        e(i+1) = y(i+1) - g(i+1);         % error signal
        K(:,i+1) = (P * V2(:,i+1)) /(lam + V2(:,i+1)' * P * V2(:,i+1));
        W(:,i+1) = W(:,i) + K(:,i+1) * e(i+1);
        P = (P - K(:,i+1) * V2(:,i+1)' * P)/lam;
    end
 fprintf(1,'The final filter coefficients are: \n');
 for i = 0 : length(W(:,i))-1
     fprintf(1,'w%d = %g\n',i,W(i+1,N));
 end

%    figure(1),stem(e); 
   
 figure,subplot(2,1,1),plot(x),title('training signal');
 subplot(2,1,2),plot(g,'r'),title('estimated training signal');
 stem(x),hold on,stem(g,'r');
 %Verification of channel equalisation
n = 1:N ; sig = sin(0.05* pi * (n-1));
ysig = filter(h,1,sig);
ysig = ysig' + normrnd(0,0.005,N,1);
outrx = filter(W(:,N)',1,ysig);
figure,plot(1:N,sig,1:N,outrx,'g');
legend('Transmitted','Received');
figure,plot(ysig),title('signal received at the receiver');