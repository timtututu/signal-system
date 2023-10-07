clc; clear;

[y, fs] = audioread('Sound_clip_2022.wav');

T = 1/ fs;        % sampling period
L = length(y);    % length of the signal
f = fs*(0:(L/2))/L;
  

t = (0: L- 1)* T;
subplot(7, 1, 1);
plot(t, y);
title('original');

Y = fft(y);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

filterOrder = 7;  % Order of filter
cutOffFreqLow = 275; % Cutoff frequency
cutOffFreqHi = 3000; % Cutoff frequency

cutOffFreqM1 = 100;     % Cutoff frequency
cutOffFreqM2 = 800;     % Cutoff frequency


% generating filter
[Lb, La]=butter(filterOrder, cutOffFreqLow/(fs/2), 'low');
[Hb, Ha]=butter(filterOrder, cutOffFreqHi/(fs/2), 'high');
% [Mb, Ma, Mc, Md] = butter(filterOrder, [cutOffFreqM1 cutOffFreqM2]/(fs/2));

% generate band-stop filter
% MP = zeros(length(Y), 1);
% for i = 5000: 20000
%     MP(i) = 1;
% end

piano = 1.7* filter(Lb, La, y);
subplot(7, 1, [2, 3]);
plot(t, piano);
title('piano');

man = bandpass(y, [cutOffFreqM1 cutOffFreqM2], fs);
subplot(7, 1, [4, 5]);
plot(t, man);
title('man');

violin = filter(Hb, Ha, y);
subplot(7, 1, [6, 7]);
plot(t, violin);
title('violin');


% sound(y, fs);          % original
% sound(piano, fs);      % piano
sound(man, fs);        % man
% sound(violin, fs);     % violin

