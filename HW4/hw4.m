clc;
clear;

%     109030012  卓鈺博
%     109030015  朱緯騰
%     109034030  凃光庭

%% Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1Q1

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

filterOrder = 7;		% Order of filter
cutOffFreqLow = 275;	% Cutoff frequency
cutOffFreqHi = 3000;	% Cutoff frequency

cutOffFreqM1 = 400;	    % Cutoff frequency
cutOffFreqM2 = 800;	    % Cutoff frequency


% generating filter
[Lb, La]=butter(filterOrder, cutOffFreqLow/(fs/2), 'low');
[Hb, Ha]=butter(filterOrder, cutOffFreqHi/(fs/2), 'high');

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


% Below are the commands for outputting the voice

% sound(y, fs);          % original
% sound((piano), fs);      % piano
% sound(man, fs);        % man
% sound(abs(violin), fs);     % violin


%% Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2Q2

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

filterOrder = 9;		% Order of filter
cutOffFreqLow = 350;	% Cutoff frequency
cutOffFreqHi = 6000;	% Cutoff frequency

cutOffFreqM1 = 2200;	    % Cutoff frequency
cutOffFreqM2 = 3200;	    % Cutoff frequency


% generating filter
[Lb, La]=butter(filterOrder, cutOffFreqLow/(fs/2), 'low');
[Hb, Ha]=butter(filterOrder, cutOffFreqHi/(fs/2), 'high');


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


%-----(a)-----%
noPiano = 1.2* man + violin;


%-----(b)-----%
omega_cello = -3000;
modulation_cello = exp(1i* omega_cello.* t.');
cello = 1.5* (modulation_cello).* violin;


omega_woman = 500;
modulation_woman = exp(1i* omega_woman.* t.');
woman = 1.5* (modulation_woman).* man;


% Below are the commands for outputting the voice

 sound(noPiano, fs);            % no piano
%sound(real(woman+ cello), fs);   % woman + cello




