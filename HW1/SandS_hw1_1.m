clear
clc

n = -15 : 1 : 15;%更改

x1n = exp(1j*(2*pi/3).*n);
figure(1);
subplot(3, 1, 1);
stem(n, x1n, 'filled', 'diamond','Color','r')%更改
title('x1(n)');

x2n = exp(1j*(3*pi/4).*n);
subplot(3, 1, 2);%更改
stem(n, x2n, 'filled', 'g')
title('x2(n)');

x3n = exp(1j*(1*pi/2).*n);
subplot(3, 1, 3);
stem(n, x3n, 'filled','*','Color','b')%更改
title('x3(n)');

sumn = x1n + x2n + x3n;
figure(2);
stem(n, sumn,'filled','x','Color','m')  %更改
title('x1(n)+x2(n)+x3(n)');


