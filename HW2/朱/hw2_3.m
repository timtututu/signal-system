clc; clear;

tmax = 0.5* 10^ -10;
tmin = -0.5* 10^ -10;

Fs = 150* 10^ 12;  % resolution (sampling freq.)
dt = 1/ Fs;

t = tmin: dt: tmax; 

N = (tmax - tmin)/ dt+ 1;
df = Fs/ N;
omega = (-N/ 2: N/ 2- 1)* df;

% square wave
T1 = 2* 10^ -12;
at = 1 .* (t >= -T1 & t <= T1) + 0 .* (t > tmin & t < tmax);


% gaussian pulse
t0 = 1* 10^ -12; % t0 = 1ps
bt = exp(-(t/ t0).^ 2);


subplot(2, 2, 1);
plot(t, at);
title('Square Wave (T1 = 2ps)');
grid on;
subplot(2, 2, 2);
plot(t, bt);
title('Gaussian Pulse (t_0 = 1 ps)');
grid on;

% perform convolution
a_conv_b = conv(at, bt);
y = a_conv_b.* dt;
x = (1: length(y)).* dt- 10* 10^ -11; % shifting 
subplot(2, 2, 3);
plot(x, y);
title('From convolution'); 
xlim([-10* 10^ -11, 10* 10^ -11])
grid on;


% perform Fourier Transform 
A = fftshift(fft(at));
B = fftshift(fft(bt));

AB = A.* B;
ift_AB = abs(ifftshift(ifft(AB)));

subplot(2, 2, 4);
plot(t, ift_AB.* dt);
xlim([-10* 10^ -11, 10* 10^ -11])
title('From I.F.T.');



grid on;

