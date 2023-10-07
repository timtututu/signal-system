clc;
clear;

step = 0.01;
t = -10: step: 10;

%causal system!!

%impulse response
ht = exp(-2*t) .* heaviside(t); 
figure(1);
subplot(3, 1, 1);
plot(t, ht)
title('h(t)');
ylim([0 2]);

%impulse function
delta = 0.5;
xt = 0 .* (t < 0 & t > delta) + 1 / delta .* (t >= 0 & t <= delta);
subplot(3, 1, 2);
plot(t, xt)
title('x(t), \Delta = 0.5');

xt_rev = fliplr(xt);

%找到time reverse後第一個有值的位置
for counter = 2 : length(t)
    if xt_rev(1, counter - 1) == 0 && xt_rev(1, counter) ~= 0
        rev_impulse_start_before_shift = counter;
        break;
    end
end

%doing convolution
%標記第一個不為零的位置
counter = 1;
while(1)
    if(ht(1, counter) ~= 0)
        have_value = counter;
        break;
    end
    counter = counter + 1;
end

%讓整張圖移到最左側
while(1)
    xt_rev = circshift(xt_rev, 1);
    if(xt_rev(1, (end - 1)) == 0  && xt_rev(1, end) == 0 && xt_rev(1, 1) ~= 0)
        break;
    end
end

% 相乘不為零時停下來
while(1)
    xt_rev = circshift(xt_rev, 1);

    if(sum(ht .* xt_rev, 'all') ~= 0)
        break;
    end
end

%找到time reverse後且shifting後第一個有值的位置
for counter = 2 : length(t)
    if xt_rev(1, counter - 1) == 0 && xt_rev(1, counter) ~= 0
        rev_impulse_start_after_shift = counter;
        break;
    end
end

middle = round(length(t) / 2); %找到中點位置
begin = middle - (rev_impulse_start_before_shift -...
    rev_impulse_start_after_shift); %找到起始作圖位置


con = zeros(1, length(t)); %產生零矩陣
%disp(length(t));

while(1)
    
    con(1, begin) = sum(xt_rev .* ht, 'all');

    if(con(1, begin) < 0.001) %因為不可能真的到0，所以要限定一個最小值
        break;
    end

    begin = begin + 1;
    xt_rev = circshift(xt_rev, 1);

end

figure(2);
subplot(2, 1, 1);
plot(t, con * step)
title('h(t) * x(t), \Delta = 0.5');
ylim([0 2]);
set(gca, 'xtick', (-10: 2: 10));
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%delta getting smaller
delta = 0.0025;
xt = 0 .* (t < 0 & t > delta) + 1 / delta .* (t >= 0 & t <= delta);
figure(1);
subplot(3, 1, 3);
plot(t, xt)
title('x(t), \Delta =', delta);

xt_rev = fliplr(xt);

%找到time reverse後第一個有值的位置
for counter = 2 : length(t)
    if xt_rev(1, counter - 1) == 0 && xt_rev(1, counter) ~= 0
        rev_impulse_start_before_shift = counter;
        break;
    end
end

%doing convolution
%標記第一個不為零的位置
counter = 1;
while(1)
    if(ht(1, counter) ~= 0)
        have_value = counter;
        break;
    end
    counter = counter + 1;
end

%讓整張圖移到最左側
while(1)
    xt_rev = circshift(xt_rev, 1);
    if(xt_rev(1, (end - 1)) == 0  && xt_rev(1, end) == 0 && xt_rev(1, 1) ~= 0)
        break;
    end
end

% 相乘不為零時停下來
while(1)
    xt_rev = circshift(xt_rev, 1);

    if(sum(ht .* xt_rev, 'all') ~= 0)
        break;
    end
end

%找到time reverse後且shifting後第一個有值的位置
for counter = 2 : length(t)
    if xt_rev(1, counter - 1) == 0 && xt_rev(1, counter) ~= 0
        rev_impulse_start_after_shift = counter;
        break;
    end
end

middle = round(length(t) / 2); %找到中點位置
begin = middle - (rev_impulse_start_before_shift -...
    rev_impulse_start_after_shift); %找到起始作圖位置


con = zeros(1, length(t)); %產生零矩陣
%disp(length(t));

while(1)
    
    con(1, begin) = sum(xt_rev .* ht, 'all');

    if(con(1, begin) < 0.001) %因為不可能真的到0，所以要限定一個最小值
        break;
    end

    begin = begin + 1;
    xt_rev = circshift(xt_rev, 1);

end

figure(2);
subplot(2, 1, 2);
plot(t, con * delta) % if delta < step
title('h(t) * x(t), \Delta =', delta);
ylim([0 2]);
set(gca, 'xtick', (-10: 2: 10));
grid on;


%below is the demonstration of sifting
disp('h(1.5) =')
disp(ht(1, 1 + (1.5 - (-10)) / step))

disp('sifting: ')

xt = circshift(xt, 1 + (1.5 / step));
figure(3);
plot(t, xt)
grid on;
set(gca, 'xtick', (-10: 0.5: 10));

answ = xt.* ht(1, 1 + (1.5 - (-10)) / step);
disp('h(t)(t - 1.5) =')
disp(answ(1,  (1.5 - (-10)) / step))





