clc; clear;

tmax = 30* 10^ -7;
tmin = -30* 10^ -7;

Fs = 0.8* 10^ 8;  % resolution (sampling freq.)
dt = 1/ Fs;

t = tmin: dt: tmax; 

N = (tmax - tmin)/ dt+ 1; %樣本數
df = Fs/ N;

dftest=1/N;

omega = (-N/ 2: N/ 2- 1)* df;


% gaussian pulse
t0 = 1* 10^ -6; % t0 = 1us
gt = exp(-(t/ t0).^ 2);

% cosine wave
pt = cos(2* pi* 10^ 7.* t);

% r(t)
% rt = gt.* pt.*dftest;
rt = gt.* pt;

rttest= rt.*2;


subplot(6, 2, 1);
plot(t, gt);
title('Gaussian Pulse, g(t)');
subplot(6, 2, 2);
plot(t, pt);
title('Cosine Wave, p(t)');


subplot(6, 2, [3, 4]);
plot(t, rt, 'r');
title('r(t) = g(t)p(t)');
% perform F.T.
R = abs(fftshift(fft(rt)));
subplot(6, 2, [5, 6]);
plot(omega, R, 'm');
title('R(j\omega)');


% perform F.T.
G = abs(fftshift(fft(gt)));
P = abs(fftshift(fft(pt)));

subplot(6, 2, 7);
plot(omega, G);
title('G(j\omega)');
subplot(6, 2, 8);
plot(omega, P);
title('P(j\omega)');

% Ry = conv(abs(G), abs(P)).* dftest;
Ry = conv(fftshift(fft(gt)), fftshift(fft(pt))).* dftest;
Rx = (1: length(Ry)).* df- 8e7;
subplot(6, 2, [9, 10]);
plot(Rx, abs(Ry), 'm');
xlim([-4e7 4e7])
title('G * P');

check = abs(Ry(241: 721)) - R;
subplot(6, 2, [11, 12]);
plot(omega, check);

