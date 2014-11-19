function [y,b,px]=spnlms(x,d,b,mu,sigma,alpha,px)
% [y,b,px]=spnlms(x,d,b,mu,sigma,alpha,px)
% NLMS algorithm. See Fig. 12.2. x,d=input and desired row
% vectors; b=LMS weight row vector; mu,sigma,alpha as in
% (12.8); px=row vector with past x values; may be "px=0".
% Inputs:
%    x      = input data row vector [x(1),x(2),...,x(N)].
%    d      = desired signal row vector, same length as x.
%    b      = row vector with L+1 adaptive FIR weights. 
%             H(z)=b(1)+b(2)*z^(-1)+...+b(L+1)*z^(-L).
%    mu     = convergence parameter in eq. 12.8; scalar.
%    sigma  = estimate of x^2 (updated); scalar.
%    alpha  = forgetting factor in eq. 12.8; scalar.
%    px     = past values of x row vector (may be "0").
% Outputs:  
%    y  = output row data vector with epsilon=d-y.
%    b  = updated adaptive weight row vector, length L+1.
%    px = updated px row vector; [x(N-1),...,x(N-L)].

N=length(x); L=length(b)-1;
if N~=length(d),
   error('SPNLMS: lengths of x and d row vectors not equal.');
end
if (mu<=0)|(mu>=1)|(sigma<=0)|(alpha<0)|(alpha>=1),
   error('SPNLMS: mu, sigma, or alpha out of range.');
end
y=zeros(1,N);
if(length(px)<L),
   px=[px,zeros(1,L-length(px))];
end
px=[0,px];
for k=1:N,
   px(1)=x(k);
   y(k)=b*px';
   if abs(y(k))>1e10,
      fprintf('\nSPNLMS warning: |y| output > 1e10.\n');
      y(k+1:N)=zeros(1,N-k);
      return
   end
   e=d(k)-y(k);
   sigma=alpha*(px(1)^2)+(1-alpha)*sigma;
   tmp=2*mu/((L+1)*sigma);
   b=b+tmp*e*px;
   px(L+1:-1:2)=px(L:-1:1);
end
px=px(2:L+1);
return
