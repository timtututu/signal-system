clear
clc

n = -10: 10;
xn = n .* (n >= -2 & n <= 2) + 0 .* (n > 2 & n < -2);
xn_rev = -xn
xn_rev = circshift(xn_rev, -5); %讓xn_rev移動到有overlap的前一刻

% figure(1);
% subplot(2, 1, 1);
% stem(n, xn)
% 
% subplot(2, 1, 2);
% stem(n, xn_rev)

counter = 9; %總共會需要做9次平移
a = 7;       %第一個stem會出現在第7個
autocon = zeros(1, 21);

while counter > 0
    xn_rev = circshift(xn_rev, 1);
    autocon(1, a) = sum(xn_rev .* xn, 'all');
    a = a + 1;
    counter = counter - 1;
end

figure(2);
%plot(, autocon)



% figure(3);
% stem(conv(xn, xn))
