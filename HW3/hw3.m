clc;
clear;

%     109030012  卓鈺博
%     109030015  朱緯騰
%     109034030  凃光庭

%% Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1

clc; clear;

Fs = 5* 1e8;%time domain resolution
delta_omega_a = 2* pi* 1e9;
delta_omega_c = 2* pi* 1e10;
L = 500;
d_omega = Fs* 2* pi;
omega = d_omega.* (-L: L);
dt = (2* pi)/ (Fs* 2* pi.* (2* L+ 1));
t = (-L: L).* dt;

%%%(a)%%%

% plotting A(jw)
A_a = zeros(1, length(omega));
for i = 1: length(omega)
    if abs(omega(i)/ delta_omega_a- 2)< 0.01
        A_a(i) = 1;
    elseif abs(omega(i)/ delta_omega_a- 1)< 0.01
        A_a(i) = 1;
    elseif abs(omega(i)/ delta_omega_a- 0)< 0.01
        A_a(i) = 1;
    elseif abs(omega(i)/ delta_omega_a+ 1)< 0.01
        A_a(i) = 1;
    elseif abs(omega(i)/ delta_omega_a+ 2)< 0.01
        A_a(i) = 1;
    end

end

figure(1);

subplot(4, 2, 1);
stem(omega./ delta_omega_a, A_a);
xlim([-2.5, 2.5]);
title('A(j\omega), \Delta\omega = 2\pi*10^9');
xlabel('\omega/\Delta\omega');

% plotting a(t)
a_a = ifftshift(ifft(ifftshift(A_a)));
subplot(4, 2, 3);
plot(t, (a_a));
title('a_1(t)');
xlabel('t');

%%%(b)%%%

% plotting A(jw)
A_b = zeros(1, length(omega));
for i = 1: length(omega)
    if abs(omega(i)/ delta_omega_a- 2)< 0.01
        A_b(i) = 1;
    elseif abs(omega(i)/ delta_omega_a- 1)< 0.01
        A_b(i) = 1;
    elseif abs(omega(i)/ delta_omega_a- 0)< 0.01
        A_b(i) = 1;
    elseif abs(omega(i)/ delta_omega_a+ 1)< 0.01
        A_b(i) = 1;
    elseif abs(omega(i)/ delta_omega_a+ 2)< 0.01
        A_b(i) = 1;
    end
end

A_b = circshift(A_b, 20);

subplot(4, 2, 2);
stem((omega./ delta_omega_a), A_b, 'r');
xlim([7.5, 12.5]);
title('A(j\omega)');
xlabel('\omega/\Delta\omega');

% plotting a(t)
a_b = ifftshift(ifft(ifftshift(A_b)));
subplot(4, 2, 4);
plot(t, (a_b), 'r');
title('a_2(t)');
xlabel('t');

%%%(c)%%%

I_a = abs(a_a).^ 2;
I_b = abs(a_b).^ 2;

subplot(4, 2, [5, 6]);
plot(t, I_a, 'g');
title('I_1(t)');
xlabel('t');
subplot(4, 2, [7, 8]);
plot(t, (I_b), 'g');
title('I_2(t)');
xlabel('t');


%% Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2

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


%% Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3Q3

clc; clear;

tmax = 15;
tmin = 0;

Fs = 8e3;  % resolution (sampling freq.)
dt = 1/ Fs;

omega_n = 100;

t = (tmin: dt: tmax)./ omega_n; 

N = (tmax - tmin)/ dt+ 1;
df = Fs/ N;
omega = (-N/ 2: N/ 2- 1)* df;

zeta1 = 0.5;
zeta2 = 1;
zeta3 = 1.5;



c1_1 = -zeta1.* omega_n + omega_n.* sqrt(zeta1.^ 2- 1);
c2_1 = -zeta1.* omega_n - omega_n.* sqrt(zeta1.^ 2- 1);
M_1 = omega_n./ (2* sqrt(zeta1.^ 2- 1));


c1_3 = -zeta3.* omega_n + omega_n.* sqrt(zeta3.^ 2- 1);
c2_3 = -zeta3.* omega_n - omega_n.* sqrt(zeta3.^ 2- 1);
M_3 = omega_n./ (2* sqrt(zeta3.^ 2- 1));


ht_1 = M_1.*(exp(c1_1.* t)- exp(c2_1.* t));
ht_2 = omega_n.^ 2.* t.* exp(-omega_n.* t);
ht_3 = M_3.*(exp(c1_3.* t)- exp(c2_3.* t));



subplot(3, 3, [1, 2, 3]);
plot(t, ht_1./ omega_n, 'r', t, ht_2./ omega_n, 'g', t, ht_3./ omega_n, 'b');
grid on;
ylabel('h(t)/\omega_n');
xlabel('t (1/\omega_n)');
legend('\zeta = 0.5', '\zeta = 1', '\zeta = 1.5');


H_1 = (M_1./ (1i.* omega- c1_1))- (M_1./ (1i.* omega- c2_1));
H_2 = (omega_n.^ 2)./ ((1i.* omega+ omega_n).^ 2);
H_3 = (M_3./ (1i.* omega- c1_3))- (M_3./ (1i.* omega- c2_3));




subplot(3, 3, [4, 5, 6]);
semilogx(omega, 20*log10(abs(H_1)), 'r', omega, 20*log10(abs(H_2)) ...
    , 'g', omega, 20*log10(abs(H_3)), 'b');
title('20log_{10}|H(j\omega)|');
xlim([1, 1e4]);
ylim([-30, 10]);
xlabel('\omega');
ylabel('dB');
legend('\zeta = 0.5', '\zeta = 1', '\zeta = 1.5');
grid on;


ph_1 = atan(imag(H_1)./ real(H_1)).* (omega <= omega_n)+ ...
    (atan(imag(H_1)./ real(H_1)) - pi).* (omega > omega_n);
ph_2 = atan(imag(H_2)./ real(H_2)).* (omega <= omega_n)+ ...
    (atan(imag(H_2)./ real(H_2)) - pi).* (omega > omega_n);
ph_3 = atan(imag(H_3)./ real(H_3)).* (omega <= omega_n)+ ...
    (atan(imag(H_3)./ real(H_3)) - pi).* (omega > omega_n);


subplot(3, 3, [7, 8, 9]);
semilogx(omega, ph_1, 'r', omega, ph_2 ...
    , 'g', omega, ph_3, 'b');
title('\angle H(j\omega)');
xlabel('\omega');
xlim([1, 1e4]);
ylim([-pi, pi/ 4]);
legend('\zeta = 0.5', '\zeta = 1', '\zeta = 1.5');
grid on;

%% Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4Q4

clc; clear;

nmax = 15;
nmin = -15;

Fs = 10;  % resolution (sampling freq.)
dn = 1/ Fs;
n = nmin: dn: nmax; 
N = (nmax - nmin)/ dn+ 1; % number of data
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

A = 1- (r.* exp(1i.* (theta2- omega)));
B = 1- (r.* exp(-1i.* (theta2+ omega)));
H1 = 1./ (1- r.* exp(-1i.* omega)).^ 2;
H2 = 1./ (A.* B);
H3 = 1./ (1+ r.* exp(-1i.* omega)).^ 2;

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


subplot(3, 3, [7, 8, 9]);
plot(omega, pH1, 'r', omega, pH2, 'g', omega, pH3, 'b');
title('\angle H(e^{j\omega})');
xlim([-pi, pi]);
xlabel('\omega', FontSize = 15);
legend('\theta = 0', '\theta = 0.5\pi', '\theta = \pi');





