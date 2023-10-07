clc; clear;

Fs = 1* 1e8;
L = 1000;
omega = Fs* 2* pi.* (-L: L);
dt = (2* pi)/ (Fs* 2* pi.* (2* L+ 1));
t = (-L: L).* dt;

beta_a = 2* pi;
beta_b = 4* pi;

omega_rep = 2* pi* 1e9;
e_pm_a = exp(1i* beta_a* sin(omega_rep.* t));
e_pm_b = exp(1i* beta_b* sin(omega_rep.* t));

E_pm_a = fftshift(fft(fftshift(e_pm_a)));
E_pm_b = fftshift(fft(fftshift(e_pm_b)));

%%%(a)%%%
figure(1);

subplot(5, 2, 1);
plot(t, ((e_pm_a)));
title('e_{PM}(t), \beta=2\pi');
xlabel('t');

subplot(5, 2, 3);
plot(omega, (abs(E_pm_a)));
title('|E_{PM}(j\omega)|');
xlabel('\omega');

subplot(5, 2, 5);
plot(omega, angle(E_pm_a));
title('\angleE_{PM}(j\omega)');
xlabel('\omega');

subplot(5, 2, 7);
plot(t, abs((e_pm_a)).^ 2);
ylim([0, 2]);
title('|e_{PM}(t)|^2');
xlabel('t');

%%%(b)%%%
subplot(5, 2, 2);
plot(t, ((e_pm_b)), 'r');
title('e_{PM}(t), \beta=4\pi');
xlabel('t');

subplot(5, 2, 4);
plot(omega, (abs(E_pm_b)), 'r');
title('|E_{PM}(j\omega)|');
xlabel('\omega');

subplot(5, 2, 6);
plot(omega, angle(E_pm_b), 'r');
title('\angleE_{PM}(j\omega)');
xlabel('\omega');

subplot(5, 2, 8);
plot(t, abs((e_pm_b)).^ 2, 'r');
ylim([0, 2]);
title('|e_{PM}(t)|^2');
xlabel('t');

%%%(c)%%%
E_pm_a_flat = abs(E_pm_a);
E_pm_b_flat = abs(E_pm_b);

e_pm_a_flat = ifftshift(ifft(ifftshift(E_pm_a_flat)));
e_pm_b_flat = ifftshift(ifft(ifftshift(E_pm_b_flat)));

subplot(5, 2, 9);
plot(t, abs((e_pm_a_flat)).^ 2);
title('|e_{PM}(t)|^2 with flat phase');
xlabel('t');

subplot(5, 2, 10);
plot(t, abs((e_pm_b_flat)).^ 2, 'r');
title('|e_{PM}(t)|^2 with flat phase');
xlabel('t');








