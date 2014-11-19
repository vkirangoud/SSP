% % /**************STATISTICAL SIGNAL PROCESSING ASSIGNMENT.*************/
%               ^^^^^^^^^^^ ^^^^^^ ^^^^^^^^^^ ^^^^^^^^^^.
%2. Channel Equalisation.
% h[n] = 0.1 * (0.5)^n, n =0,1,...8.
%x training signal.
clear;
clc;
N = 2000;
% Generate bit stream

x = rand(N,1);
z0 = find(x < 0.5);
z1 = find(x >= 0.5);
x(z0) = -1*ones(size(z0));
x(z1) = +1*ones(size(z1));
    
n=1:8; h = 0.1 * ((0.5) .^ n);
y = filter(h,1,x);
% yn = y + sqrt(0.005) * randn(1,N);
yn = y + normrnd(0,0.005,N,1);
L = 20;
mu = 0.25;
% % calculation of range of step size mu .
% % ----------- -- ----- -- ---- ---- --
% c = yn(1:L)';
% rv2 = 1/N * toeplitz(c,yn) * yn;
% evalues = eig(toeplitz(rv2));
% mumax = 2/max(evalues);
% fprintf(1,'\n\nenter the value of mu in the range 0 < mu < %g : \n', mumax);
% mu = input('mu = ');
% %---------------------------------------------------------------------
V2 = flipud(buffer(yn,L,L-1)); % A buffer whose columns will act as data for each iteration. 
W(:,1) = zeros(L,1);
    for i = 1 : N
          g(i) = W(:,i)' * V2(:,i);    % output at time i.
        e(i) = x(i) - g(i);         % error signal
 
        W(:,i+1) = W(:,i) + mu * e(i) * V2(:,i) / (norm(V2(:,i),2)^2);
    end
  fprintf(1,'The final filter coefficients are: \n');
 for i = 0 : length(W(:,i))-1
     fprintf(1,'w%d = %g\n',i,W(i+1,N));
 end

          
   
 subplot(2,1,1),plot(x),title('training signal');
 subplot(2,1,2),plot(g,'r'),title('estimated training signal');