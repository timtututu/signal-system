clc;
clear;

%     109030012  卓鈺博
%     109030015  朱緯騰
%     109034030  凃光庭

%% Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1

clc; clear;


tmax = 10* 10^ -10;
tmin = -10* 10^ -10;

Fs = 5* 10^ 12;  % resolution (sampling freq.)

t = tmin: 1/ Fs: tmax; 

N = (tmax - tmin)/ (1/ Fs)+ 1; % the number of sampling in freq. domain
df = Fs/ N;
omega = (-N/ 2: N/ 2- 1)* df* 2* pi;

t0_1 = 10* 10^ -12; % t0 = 10ps
t0_2 = 1* 10^ -12; % t0 = 1ps

% generate Gaussian pulse
gp_1 = exp(-(t/ t0_1).^ 2);
gp_2 = exp(-(t/ t0_2).^ 2);

figure(1);
subplot(2, 2, 1);
plot(t, gp_1);
ylim([0, 1.2]);
title('t_0 = 10 ps');
xlabel('t(s)');
subplot(2, 2, 2);
plot(t, gp_2);
ylim([0, 1.2]);
title('t_0 = 1 ps');
xlabel('t(s)');

fft_1 = fftshift(fft(fftshift(gp_1)));
fft_2 = fftshift(fft(fftshift(gp_2)));

subplot(2, 2, 3);
plot(omega, fft_1);
xlim([-0.05e14, 0.05e14]);
ylim([0, 100]);
title('After FFT');
xlabel('\omega(s^{-1})');
subplot(2, 2, 4);
plot(omega, fft_2);
ylim([0, 10]);
title('After FFT');
xlabel('\omega(s^{-1})');

%% Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2

clc; clear;

tmax = 15* 10^ -10;
tmin = -15* 10^ -10;

Fs = 7* 10^ 12;  % resolution (sampling freq.)

t = tmin: 1/ Fs: tmax; 

N = (tmax - tmin)/ (1/ Fs)+ 1; % the number of sampling in freq. domain
df = Fs/ N;
omega = (-N/ 2: N/ 2- 1)* df* 2* pi;


%______________the signal in 1st question________________

t0_1 = 10* 10^ -12; % t0 = 10ps
t0_2 = 1* 10^ -12; % t0 = 1ps

gp_1 = exp(-(t/ t0_1).^ 2);
gp_2 = exp(-(t/ t0_2).^ 2);

figure(1);
subplot(2, 2, 1);
plot(t, gp_1);
title('t_0 = 10 ps');
ylim([0, 2]);
xlim([-3e-11, 3e-11]);
xlabel('t(s)');
grid on;
subplot(2, 2, 2);
plot(t, gp_2);
title('t_0 = 1 ps');
ylim([0, 2]);
xlim([-3e-11, 3e-11]);
xlabel('t(s)');
grid on;

%______________the signal in 1st question________________

t0 = 2e-12;  % time shift
% linear spectral phase
lsp = exp(1i.* omega.* t0);

fft_1_lsp = fftshift(lsp.* fft(fftshift(gp_1)));
fft_2_lsp = fftshift(lsp.* fft(fftshift(gp_2)));

ifft_1 = ifftshift(ifft(ifftshift(fft_1_lsp)));
ifft_2 = ifftshift(ifft(ifftshift(fft_2_lsp)));

subplot(2, 2, 3);
plot(t, ifft_1);
title('After Shifting');
ylim([0, 2]);
xlim([-3e-11, 3e-11]);
xlabel('t(s)');
xline((-2* 10^ -12), '-', 'left shift 2 ps');
grid on;
subplot(2, 2, 4);
plot(t, ifft_2);
title('After Shifting');
ylim([0, 2]);
xlim([-3e-11, 3e-11]);
xlabel('t(s)');
xline((-2e-12), '-', 'left shift 2 ps');
grid on;

%% Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3

clc; clear;

tmax = 0.5* 10^ -10;
tmin = -0.5* 10^ -10;

Fs = 25* 10^ 12;  % resolution (sampling freq.)
dt = 1/ Fs;

t = tmin: 1/ Fs: tmax; 

% square wave
T1 = 2* 10^ -12;
at = 1 .* (t >= -T1 & t <= T1) + 0 .* (t > tmin & t < tmax);


% gaussian pulse
t0 = 1* 10^ -12; % t0 = 1ps
bt = exp(-(t/ t0).^ 2);


subplot(3, 2, 1);
plot(t, at);
title('a(t), Square Wave (T1 = 2ps)');
xlabel('t(s)');
grid on;
subplot(3, 2, 2);
plot(t, bt);
title('b(t), Gaussian Pulse (t_0 = 1 ps)');
xlabel('t(s)');
grid on;

% perform convolution
a_conv_b = conv(at, bt);
y = a_conv_b.* dt;
x = (1: length(y)).* dt- 10* 10^ -11; % shifting 
subplot(3, 2, 3);
plot(x, y, 'r');
title('c(t) = a(t) * b(t)'); 
xlabel('t(s)');
xlim([-5* 10^ -11, 5* 10^ -11])
grid on;


% perform Fourier Transform 
A = fftshift(fft(fftshift(at)));
B = fftshift(fft(fftshift(bt)));

AB = A.* B;
ift_AB = ifftshift(ifft(ifftshift(AB)));

subplot(3, 2, 4);
plot(t, ift_AB.* dt, 'r');
title('c(t) = F^{-1}\{AB\}');
xlabel('t(s)');
grid on;

subplot(3, 2, [5, 6]);
check = y(length(y)/ 4+ 2: 3* length(y)/ 4+ 2) - ift_AB.* dt;
plot(t, check, 'g');
title('a(t) * b(t) - F^{-1}\{AB\}, to prove that two results are the same.');
xlabel('t(s)');
grid on;

%% Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4

clc; clear;

tmax = 900* 10^ -7;
tmin = -900* 10^ -7;

Fs = 0.8* 10^ 8;  % resolution (sampling freq.)
dt = 1/ Fs;

t = tmin: dt: tmax; 

N = (tmax - tmin)/ dt+ 1;
df = Fs/ N;
omega = (-N/ 2: N/ 2- 1)* df;


% gaussian pulse
t0 = 1* 10^ -6; % t0 = 1us
gt = exp(-(t/ t0).^ 2);

% cosine wave
pt = cos(2* pi* 10^ 7.* t);

% r(t)
rt = gt.* pt;


subplot(5, 2, 1);
plot(t, gt);
title('Gaussian Pulse, g(t)');
xlabel('t(s)');
subplot(5, 2, 2);
plot(t, pt);
title('Cosine Wave, p(t)');
xlabel('t(s)');


subplot(5, 2, [3, 4]);
plot(t, rt, 'r');
title('r(t) = g(t)p(t)');
xlabel('t(s)');
% perform F.T.
R = fftshift(fft(fftshift(rt)));
subplot(5, 2, [5, 6]);
plot(omega, R, 'm');
xlim([-3e7 3e7]);
title('R(j\omega)');
xlabel('\omega(s^{-1})');


% perform F.T.
G = fftshift(fft(fftshift(gt)));
P = fftshift(fft(fftshift(pt)));


Ry = conv(G, P)/ N;
Rx = (1: length(Ry)).* df- 8e7;
subplot(5, 2, [7, 8]);
plot(Rx, Ry, 'm');
xlim([-3e7 3e7]); 
title('(1/2\pi) G(j\omega) * P(j\omega)');
xlabel('\omega(s^{-1})');


check = Ry(length(Ry)/ 4+ 1: 3* length(Ry)/ 4+ 1) - R;
subplot(5, 2, [9, 10]);
plot(omega, check, 'g');
xlim([-3e7 3e7]);
title(['(1/2\pi) G(j\omega) * P(j\omega) - R(j\omega), ' ...
    'to prove that two results are the same.']);
xlabel('\omega(s^{-1})');







































