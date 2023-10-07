clc; clear;

tmax = 15;
tmin = 0;

Fs = 1e3;  % resolution (sampling freq.)
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







