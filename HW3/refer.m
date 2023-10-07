clc;clear all;
%%
%Q1 
%Gaussian Pulse t0 values of 1 ps
t0 = 1e-12;
dt = t0/10;
L = 100;
t = (-L:L)*dt;
Fs = 1/((2*L+1)*dt);

figure(111)
g =exp(-((t./t0).^2));
plot(t,g)
hold on
title('Gaussian Pulse t0 values of 10 ps')
xlabel('t(s)')
ylabel('g(t)')
grid on
n = length(t);
g_FFT = fft(fftshift(g),n);
g_FFTShift = fftshift(g_FFT);
w = 2*pi*(-L:L)*Fs;
g_jw = t0*(pi^0.5)*exp(-((w*t0).^2)/4);
figure(112)
plot(w,g_FFTShift*dt,'LineWidth',2)
figure(112)
hold on
plot(w,g_jw,'o')
xlabel('w')
ylabel('G(jw)')
grid on
%Gaussian Pulse t0 values of 10 ps
t0 = 10e-12;
dt = t0/10;
L = 100;
t = (-L:L)*dt;
Fs = 1/((2*L+1)*dt);
%第二章
figure(121)
g =exp(-((t./t0).^2));
plot(t,g)
hold on
title('Gaussian Pulse t0 values of 10 ps')
xlabel('t(s)')
ylabel('g(t)')
grid on
n = length(t);
g_FFT = fft(fftshift(g),n);
g_FFTShift = fftshift(g_FFT);
w = 2*pi*(-L:L)*Fs;
g_jw = t0*(pi^0.5)*exp(-((w*t0).^2)/4);
figure(122)
plot(w,g_FFTShift*dt,'LineWidth',2)
figure(122)
hold on
plot(w,g_jw,'o')
xlabel('w')
ylabel('G(jw)')
grid on

%%
%Q2
t0=1e-12;
Fs = 1e14; % 採樣頻率
dt = 1/Fs; % 採樣周期
L = 1000; % 信號長度
t = (-L:L)*dt; % Time span
t_shift = -2e-12;

n = length(t);
w = 2*pi*(-L:L)*Fs/n;
g_jw = exp(-1i.*w.*t_shift).*t0*(pi^0.5).*exp(-((w.*t0).^2)/4);

g_FFT = ifft(ifftshift(g_jw),n);
g_FFTShift = ifftshift(g_FFT);

figure(21)
plot(t,g_FFTShift,'LineWidth',2)
grid on