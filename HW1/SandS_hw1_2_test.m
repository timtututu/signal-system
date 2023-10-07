clear
clc

%xn = [-2: 2];
%figure(1);
%grid on;
%stem(xn, xn)

%figure(2);
%grid on;
%stem(conv(xn, xn))
%xn= zeros(1,3);
%xn = [-2: 2]
%xn(1:2)=[2:-1:1]
n = -5: 5;
xn = n .* (n >= -2 & n <= 2) + 0 .* (n > 2 & n < -2);
xn_rev = fliplr(xn)