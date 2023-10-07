clc; clear;

nmax = 15;
nmin = -15;

Fs = 10;  % resolution (sampling freq.)
dn = 1/ Fs;
n = nmin: dn: nmax;
N = (nmax - nmin)/ dn+ 1;
df = Fs/ N;
omega = (-N/ 2: N/ 2- 1)* df;

Fs_1 = 1;  % resolution (sampling freq.)
dn_1 = 1/ Fs_1;
n_1 = nmin: dn_1: nmax; 

r = 0.5;
theta1 = 0;
theta2 = 0.5* pi;
theta3 = 1* pi;

% Time Domain
h1 = 0.* (n< 0)+ ((n+ 1).* (r.^ n)).* (n>= 0);
h2 = 0.* (n< 0)+ ...
    ((r.^ n).* ((sin((n+ 1).* theta2))./ (sin(theta2)))).* (n>= 0);
h3 = 0.* (n< 0)+ ((n+ 1).* ((-r).^ n)).* (n>= 0);

h1_1 = 0.* (n_1< 0)+ ((n_1+ 1).* (r.^ n_1)).* (n_1>= 0);
h2_1 = 0.* (n_1< 0)+ ...
    ((r.^ n_1).* ((sin((n_1+ 1).* theta2))./ (sin(theta2)))).* (n_1>= 0);
h3_1 = 0.* (n_1< 0)+ ((n_1+ 1).* ((-r).^ n_1)).* (n_1>= 0);

figure(1);

subplot(3, 3, 1);
stem(n_1, h1_1, 'r');
title('r = 0.5, \theta = 0', FontSize = 15);
hold on;
subplot(3, 3, 2);
stem(n_1, h2_1, 'g');
title('r = 0.5, \theta = 0.5\pi', FontSize = 15);
hold on;
subplot(3, 3, 3);
stem(n_1, h3_1, 'b');
title('r = 0.5, \theta = \pi', FontSize = 15);

% subplot(3, 3, [1, 2, 3]);
% plot(n_1, h1_1, 'r', n_1, h2_1, 'g', n_1, h3_1, 'b');
% title('\angle H(e^{j\omega})');
% xlim([-pi, pi]);
% xlabel('\omega', FontSize = 15);
% legend('\theta = 0', '\theta = 0.5\pi', '\theta = \pi');


% figure(2);
% Frequency Domain


% A = (exp(1i* theta2))/ (2i* sin(theta2));
% B = (exp(-1i* theta2))/ (2i* sin(theta2));
% H2 = (A./ (1- r* exp(1i* theta2).* exp(-1i.* omega))) ...
%     + (B./ (1- r* exp(-1i* theta2).* exp(-1i.* omega)));

A = 1- (r.* exp(1i.* (theta2- omega)));
B = 1- (r.* exp(-1i.* (theta2+ omega)));
H1 = 1./ (1- r.* exp(-1i.* omega)).^ 2;
H2 = 1./ (A.* B);
H3 = 1./ (1+ r.* exp(-1i.* omega)).^ 2;

% H2 = fftshift(fft(fftshift(h2)));

% subplot(3, 3, 4);
% plot(omega, 20.* log10(abs(H1)));
% title('20log_{10}|H(e^{j\omega})|');
% xlabel('\omega', FontSize = 15);
% ylabel('dB');
% xlim([-pi, pi]);
% ylim([-15, 15]);
% subplot(3, 3, 5);
% plot(omega, 20.* log10(abs(H2)));
% title('20log_{10}|H(e^{j\omega})|');
% xlabel('\omega', FontSize = 15);
% ylabel('dB');
% xlim([-pi, pi]);
% ylim([-15, 15]);
% subplot(3, 3, 6);
% plot(omega, 20.* log10(abs(H3)));
% title('20log_{10}|H(e^{j\omega})|');
% xlabel('\omega', FontSize = 15);
% ylabel('dB');
% xlim([-pi, pi]);
% ylim([-15, 15]);


subplot(3, 3, [4, 5, 6]);
plot(omega, 20.* log10(abs(H1)), 'r', omega, 20.* log10(abs(H2)), 'g', ...
    omega, 20.* log10(abs(H3)), 'b');
title('20log_{10}|H(e^{j\omega})|');
xlabel('\omega', FontSize = 15);
ylabel('dB');
xlim([-pi, pi]);
legend('\theta = 0', '\theta = 0.5\pi', '\theta = \pi');


pH1 = atan(imag(H1)./ real(H1));
pH2 = atan(imag(H2)./ real(H2));
pH3 = atan(imag(H3)./ real(H3));


% subplot(3, 3, 7);
% plot(omega, pH1);
% title('\angle H(e^{j\omega})');
% xlabel('\omega', FontSize = 15);
% xlim([-pi, pi]);
% ylim([-pi/ 2, pi/ 2]);
% 
% subplot(3, 3, 8);
% plot(omega, pH2);
% title('\angle H(e^{j\omega})');
% xlabel('\omega', FontSize = 15);
% xlim([-pi, pi]);
% ylim([-pi/ 2, pi/ 2]);
% 
% subplot(3, 3, 9);
% plot(omega, pH3);
% title('\angle H(e^{j\omega})');
% xlabel('\omega', FontSize = 15);
% xlim([-pi, pi]);
% ylim([-pi/ 2, pi/ 2]);

subplot(3, 3, [7, 8, 9]);
plot(omega, pH1, 'r', omega, pH2, 'g', omega, pH3, 'b');
title('\angle H(e^{j\omega})');
xlim([-pi, pi]);
xlabel('\omega', FontSize = 15);
legend('\theta = 0', '\theta = 0.5\pi', '\theta = \pi');



















