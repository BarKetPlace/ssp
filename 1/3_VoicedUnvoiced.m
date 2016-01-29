%% 3- Voiced and Unvoiced Speech Sounds
clear all
close all
clc

load assignment1.mat
S = male_short;
Fs = 8000;
scaling_f = 1;
mute = 0;

figure, subplot(2,1,1)
spectrogram(S, 200, 50, 200, Fs, 'yaxis'); colorbar('off');
title(['male short Frequency vs time']);
subplot(2,1,2);
display_(S, Fs, scaling_f, mute);
%%

%% DFT
close all

% x is a vector of speech, N is the frame length (must be even),
% S is the first sample of the frame to be analyzed

%UnVoiced : [f] [k] [p][s][t][tch][th][ss]
% Sound s from the word "moSt"

Suv=5900;
x = male_short;
Nuv = 1000;
xfuv = x(Suvcd:Suv+Nuv-1).*hanning(Nuv);
Xuv = fft(xfuv);
figure,
plot(10*log10(abs(Xuv(1:Nuv/2+1)).^2));
title('UnVoiced');
%Voiced : [v][g][b][z][d][>th<is]

Sv=122;
Nv =450;
xfv = x(Sv:Sv+Nv-1).*hanning(Nv);
Xv = fft(xfv);
figure, 
plot(10*log10(abs(Xv(1:Nv/2+1)).^2));
title('Voiced');