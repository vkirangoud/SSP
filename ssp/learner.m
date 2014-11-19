echo on
% LearnER.m - Example of Error Correcting Learning
% copyright (c) 1996-2000 by Yu Hen Hu
% Created: 9/2/96
% Modified: 1/28/2000
% Modified: 9/3/2001 add additional runs of LMS and display weight converge curve
% 
clear all
x=[  1     1    1    1
    0.5 -0.4  1.1  0.7
    0.8  0.4 -0.3  1.2];
d=[1    0    0    1];
w=zeros(3,1);
weight=[];
eta=.01;
echo off; pause
for n=1:4,
   y(n) = w'*x(:,n);
   e(n) = d(n) - y(n);
   w=w+eta*e(n)*x(:,n);
   weight=[weight w];
   ['iteration #' int2str(n) ':']
   weight
   pause
end
for m=1:499,
   x0=randomize(x')'; % change the order of presentation of x
   for n=1:4,
      y(n) = w'*x0(:,n);
      e(n) = d(n) - y(n);
      w=w+eta*e(n)*x0(:,n);
      weight=[weight w];
   end
end
figure(1),
subplot(311),plot([1:size(weight,2)],weight(1,:)),ylabel('w0')
title('convergence curve of the weights')
subplot(312),plot([1:size(weight,2)],weight(2,:)),ylabel('w1')
subplot(313),plot([1:size(weight,2)],weight(3,:)),ylabel('w2')

echo on
% Batched mode LS solution
R = x*x'
rho=sum((x*diag(d))')'
w_ls = inv(R)*rho

error = d - w_ls'*x;
ernorm = error*error'
echo off