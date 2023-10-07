clc; clear;

% resolution
t = -0.4* 10^ -10: 7* 10^ -13: 0.4* 10^ -10; % smaller
%t = -0.4* 10^ -10: 5* 10^ -15: 0.4* 10^ -10; % bigger

t0_1 = 10* 10^ -12; % t0 = 10ps
t0_2 = 1* 10^ -12; % t0 = 1ps

% generate Gaussian pulse
gp_1 = exp(-(t/ t0_1).^ 2);
gp_2 = exp(-(t/ t0_2).^ 2);

figure(1);
subplot(3, 2, 1);
plot(t, gp_1);
title('t_0 = 10 ps');
xlabel('t');
subplot(3, 2, 2);
plot(t, gp_2);
title('t_0 = 1 ps');
xlabel('t');

fft_1 = fft(gp_1);
fft_2 = fftshift(fft(gp_2));

subplot(3, 2, 3);
plot3(t, imag(fft_1), real(fft_1), ...
    t, imag(fft_1), -40.* ones(size(t)), 'r', ...
    t, ones(size(t)), real(fft_1), 'g');
title('FFT');
xlabel('\omega (freqency)');
ylabel('imaginary part');
zlabel('real part');
grid on;

subplot(3, 2, 4);
plot3(t, imag(fft_2), real(fft_2), ...
    t, imag(fft_2), -5.* ones(size(t)), 'r', ...
    t, ones(size(t)), real(fft_2), 'g');
title('FFT');
xlabel('\omega (freqency)');
ylabel('imaginary part');
zlabel('real part');
grid on;

subplot(3, 2, 5);
plot(t, fftshift_1);
title('FFTshift');
xlabel('\omega');
subplot(3, 2, 6);
plot(t, fftshift_2);
title('FFTshift');
xlabel('\omega');

% THE DIFFERENCE BETWEEN FFT / FFTSHIFT ??


