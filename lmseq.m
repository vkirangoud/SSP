% Adaptive equalizer based on LMS algorithm.
% 
% Rich Kozick, ELEC 470, Spring 1998

T = 1;      % Bit period
tau = 3;    % Time constant of channel
SNR = 100;   % Ratio of signal power to noise power (NOT in dB)

dt = 0.01;  % Sampling time in simulation
N = 250;    % Number of training bits to generate
Ndata = 100; % Number of data bits to generate

clear t1 t2 c x y

% Create output pulse: rectangular pulse convolved with first-order
% low-pass filter impulse response.

t1 = (dt:dt:T)';
c(1:100,1) = 1 - exp(-t1/tau);
t2 = (T+dt:dt:T+5*tau)';
c(101:100+length(t2),1) = c(100)*exp(-(t2-T)/tau);

% The following channel is different from RC LPF

% c(101:100+length(t2),1) = c(100)*exp(-(t2-T)/tau).*(1+((t2-T-dt).^2)/20);

figure(1)
plot([t1; t2], c)
xlabel('Time (sec)')
ylabel('c(t)')
title('Smeared pulse c(t)')

% Generate bit stream

b = rand(N,1);
z0 = find(b < 0.5);
z1 = find(b >= 0.5);
b(z0) = -1*ones(size(z0));
b(z1) = +1*ones(size(z1));

% Create received signal with ISI

nT = T/dt;
nc = length(c);
nx = N*nT;
x = zeros(nx, 1);
for n=1:N
  i1 = (n-1)*nT;
  y = [zeros(i1,1); b(n)*c; zeros(N*nT-i1-nc,1)];
  x = x + y(1:nx);
end

% Add noise of the specified level

sp = sum(x.*x)/length(x);  % Signal power
np = sp/SNR;               % Noise power
noise = sqrt(np) * randn(length(x),1);
x = x + noise;

% Plot eye diagram

figure(2)
t3 = dt:dt:2;
plot(t3, x(1:200));
hold on
for n=3:2:N
  plot(t3, x((n-1)*nT+1:(n+1)*nT));
end
hold off
xlabel('Time (sec)')
title('EYE DIAGRAM WITH NO EQUALIZATION')

% Define length of equalizer

Ne = 5*tau/T;

% Get samples of c(t) to solve for ZF equalizer weights

cT = c(nT:nT:nc);
csamp = [zeros(2*Ne,1); cT; zeros(2*Ne+1-length(cT),1)];

% Construct the matrix on the left side of ZF equalizer equation

C = zeros(2*Ne+1,2*Ne+1);
for ne = 1:2*Ne+1
  C(ne,:) = csamp(2*Ne+ne:-1:ne)';
end

% Right side of ZF equalizer weight equation

r = [zeros(Ne,1); 1; zeros(Ne,1)];

% Solve for ZF equalizer weights

wzf = C \ r;

% LMS algorithm to estimate the equalizer weights in real-time
% using the training data

mu = 0.2;
xT = x(nT:nT:nx);
w = zeros(2*Ne+1,1);     % Initialize weights to zero
wold = w;
wrec = w;
for k=(Ne+1):(N-Ne)
  xk = xT((k+Ne):-1:(k-Ne));
  yk = w'*xk;
  ek = b(k) - yk;
  w = wold + mu*ek*xk;
  wrec = [wrec, w];
  wold = w;
end

% Process the received signal with the LMS equalizer

nw = length(w);
z99 = [1; zeros(nT-1,1)];
hzf = kron(wzf, z99);        % Impulse response of ZF equalizer
yall = conv(x, hzf);         % Do the equalization filtering
yzf = yall((Ne*nT+1):(length(yall)-(Ne+1)*nT)+1);

h = kron(w, z99);        % Impulse response of LMS equalizer
yall = conv(x, h);       % Do the equalization filtering
y = yall((Ne*nT+1):(length(yall)-(Ne+1)*nT)+1);

% Eye diagrams of equalized signal

figure(3)
plot(t3, yzf(1:200));
hold on
for n=3:2:N
  plot(t3, yzf((n-1)*nT+1:(n+1)*nT));
end
hold off
xlabel('Time (sec)')
title('EYE DIAGRAM AFTER ZF EQUALIZER')

figure(4)
plot(t3, y(1:200));
hold on
for n=3:2:N
  plot(t3, y((n-1)*nT+1:(n+1)*nT));
end
hold off
xlabel('Time (sec)')
title('EYE DIAGRAM AFTER LMS EQUALIZER')

figure(5)
plot(wrec')
xlabel('TRAINING SAMPLE NUMBER')
ylabel('WEIGHT VALUE')
title('EVOLUTION OF LMS EQUALIZER WEIGHTS')

% Now that training is over, do data transmission.

% Generate bit stream

b = rand(Ndata,1);
z0 = find(b < 0.5);
z1 = find(b >= 0.5);
b(z0) = -1*ones(size(z0));
b(z1) = +1*ones(size(z1));

% Create received signal with ISI

nx = Ndata*nT;
x = zeros(nx, 1);
for n=1:Ndata
  i1 = (n-1)*nT;
  y = [zeros(i1,1); b(n)*c; zeros(Ndata*nT-i1-nc,1)];
  x = x + y(1:nx);
end

% Add noise of the specified level

noise = sqrt(np) * randn(length(x),1);
x = x + noise;

% Perform equalization

yall = conv(x, hzf);    % ZF
yzf = yall((Ne*nT+1):(length(yall)-(Ne+1)*nT)+1);

yall = conv(x, h);      % Wiener
y = yall((Ne*nT+1):(length(yall)-(Ne+1)*nT)+1);

% Compute number of bit errors

xT = x(nT:nT:nx);
dz0 = find(xT < 0);
dz1 = find(xT >= 0);

db = b;
db(dz0) = -1*ones(size(dz0));
db(dz1) = +1*ones(size(dz1));

err = find(db ~= b);

fprintf('No equalizer: %d bits out of %d in error\n', length(err), N);

yT = yzf(nT:nT:nx);
dz0 = find(yT < 0);
dz1 = find(yT >= 0);

db = b;
db(dz0) = -1*ones(size(dz0));
db(dz1) = +1*ones(size(dz1));
err = find(db ~= b);
fprintf('ZF equalizer: %d bits out of %d in error\n', length(err), Ndata);

yT = y(nT:nT:nx);
dz0 = find(yT < 0);
dz1 = find(yT >= 0);

db = b;
db(dz0) = -1*ones(size(dz0));
db(dz1) = +1*ones(size(dz1));

err = find(db ~= b);

fprintf('LMS equalizer: %d bits out of %d in error', length(err), Ndata);
fprintf(' (with %d training bits)\n', N);










