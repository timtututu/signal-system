clear
clc

n = -10: 0.01: 10;
o = 4;
p = 0;
q = -1;
r = 2;
yn = n.*(n >= p & n <= o) + 0.*(n < p & n > o);
xn = n.^3.*(n >= q & n <= r) + 0.*(n < q & n > r);
subplot(4,1,1);
plot(n, xn)

yn_rev = flip(yn);
subplot(4,1,2);
plot(n, yn_rev)

a = 0;
z = zeros(1, 2001);
while(1)
    b = -r - (o-(p))/2 + a;
    z(1, 10 + b) = sum(xn.*circshift(yn_rev, b), 'all');
    a = a + 1;
    if(a >= 1995)
        break;
    end
end
subplot(4,1,3);
plot(n, z * 1)

subplot(4,1,4);
plot(conv(xn, yn * 1))
