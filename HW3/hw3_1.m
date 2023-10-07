clear;

fmax = 200* pi* 1e9;
fmin = -200* pi* 1e9;
df = 0.4* pi* 1e9;

delta_omega_a = 2* pi* 1e9;
delta_omega_b = 2* pi* 1e10;

omega = (fmin: df: fmax);

%%%(a)%%%

% plotting A(jw)
int_index_0_a = abs(mod(omega./ delta_omega_a, 1) - 0) <= 0.01;
int_index_1_a = abs(mod(omega./ delta_omega_a, 1) - 1) <= 0.01;
A_a = zeros(1, length(omega));
A_a(int_index_0_a) = 1;
A_a(int_index_1_a) = 1;

figure(1);

subplot(4, 2, 1);
stem(omega./ delta_omega_a, A_a);
xlim([-10, 10]);
title('A(j\omega), \Delta\omega = 2\pi*10^9');
xlabel('\omega/\Delta\omega');

% plotting a(t)
a_a = ifftshift(ifft(ifftshift(A_a)));
subplot(4, 2, 3);
plot(omega, abs(a_a));
title('a_1(t)');
xlabel('t');

% remember to revise the time scale

%%%(b)%%%

% plotting A(jw)
int_index_0_b = abs(mod(omega./ delta_omega_b, 1) - 0) <= 0.01;
int_index_1_b = abs(mod(omega./ delta_omega_b, 1) - 1) <= 0.01;
A_b = zeros(1, length(omega));
A_b(int_index_0_b) = 1;
A_b(int_index_1_b) = 1;

subplot(4, 2, 2);
stem(omega./ delta_omega_b, A_b, 'r');
xlim([-10, 10]);
title('A(j\omega), \Delta\omega = 2\pi*10^{10}');
xlabel('\omega/\Delta\omega');

% plotting a(t)
a_b = ifftshift(ifft(ifftshift(A_b)));
subplot(4, 2, 4);
plot(omega, abs(a_b), 'r');
title('a_2(t)');
xlabel('t');

%%%(c)%%%

I_a = abs(a_a).^ 2;
I_b = abs(a_b).^ 2;

subplot(4, 2, [5, 6]);
plot(omega, I_a, 'g');
title('I_1(t)');
subplot(4, 2, [7, 8]);
plot(omega, I_b, 'g');
title('I_2(t)');






