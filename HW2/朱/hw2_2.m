clc; clear;

tmax = 0.4* 10^ -10;
clc; clear;

tmax = 0.4* 10^ -10;
tmin = -0.4* 10^ -10;

Fs = 3* 10^ 12;  % resolution

t = tmin: 1/ Fs: tmax; 

N = (tmax - tmin)/ (1/ Fs)+ 1;
df = Fs/ N;
omega = 2*pi*(-N/ 2: N/ 2- 1)* df;




%______________the signal in 1st question________________

t0_1 = 10* 10^ -12; % t0 = 10ps
t0_2 = 1* 10^ -12; % t0 = 1ps

% generate Gaussian pulse
% gp_1 = exp(-((t+ 2* 10^ -12)/ t0_1).^ 2);
gp_1 = exp(-(t/ t0_1).^ 2);
gp_2 = exp(-(t/ t0_2).^ 2);

figure(1);
subplot(3, 2, 1);

plot(t, gp_1);
yline(1/2, '-', 'delta t');
title('t_0 = 10 ps');
ylim([0, 5]);
xlabel('t');
grid on;
subplot(3, 2, 2);
plot(t, gp_2);
title('t_0 = 1 ps');
ylim([0, 5]);
xlabel('t');
grid on;

fft_1 = abs(fftshift(fft(gp_1)));
fft_2 = abs(fftshift(fft(gp_2)));

% subplot(4, 2, 3);
% plot3(t, imag(fft_1), real(fft_1), ...
%     t, imag(fft_1), -50.* ones(size(t)), 'r', ...
%     t, ones(size(t)), real(fft_1), 'g');
% title('FFT');
% xlabel('\omega (freqency)');
% ylabel('imaginary part');
% zlabel('real part');
% grid on;
% 
% subplot(4, 2, 4);
% plot3(t, imag(fft_2), real(fft_2), ...
%     t, imag(fft_2), -5.* ones(size(t)), 'r', ...
%     t, ones(size(t)), real(fft_2), 'g');
% title('FFT');
% xlabel('\omega (freqency)');
% ylabel('imaginary part');
% zlabel('real part');
% grid on;

%______________the signal in 1st question________________


% omega = -0.4* 10^ 1: 3* 10^ 1: 0.4* 10^ 1;

% omega = 1./ t;

t0 = 2e-12;  % time shift
% linear spectral phase
% lsp = omega;
lsp = exp(1i.* omega.* t0);

fft_1_lsp = lsp.* fftshift(fft(gp_1));
fft_2_lsp = lsp.* fftshift(fft(gp_2));

subplot(3, 2, 3);
plot(omega, abs(fft_1_lsp), 'r');
yline(20, '-', 'BW');
xlabel('\omega (freqency)');
title('FFT (after linear spectral phase)');
grid on;
subplot(3, 2, 4);
plot(omega, abs(fft_2_lsp), 'r');
xlabel('\omega (freqency)');
title('FFT (after linear spectral phase)');
grid on;

ifft_1 = abs(ifft(fft_1_lsp));
ifft_2 = abs(ifft(fft_2_lsp));

subplot(3, 2, 5);
plot(t, ifft_1);
title('After IFFT');
ylim([0, 3]);
xlabel('t');
xline((-2* 10^ -12), '-', 'left shift 2 ps');
grid on;
subplot(3, 2, 6);
plot(t, ifft_2);
title('After IFFT');
ylim([0, 3]);
xlabel('t');
xline((-2e-12), '-', 'left shift 2 ps');
grid on;

% % To prove the same
% prove_1 = gp_1 - ifft_1;
% prove_2 = gp_2 - ifft_2;
% 
% subplot(4, 2, 7);
% plot(t, prove_1);
% % title('After IFFT');
% ylim([0, 1]);
% xlabel('t');
% subplot(4, 2, 8);
% plot(t, prove_2);
% % title('After IFFT');
% ylim([0, 1]);
% xlabel('t');




% Why abs()?


