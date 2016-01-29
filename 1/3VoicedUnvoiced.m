%% 3- Voiced and Unvoiced Speech Sounds
clear all
close all
clc

load assignment1.mat
S = male_short;
Fs = 8000;
scaling_f = Fs;
mute = 0;
% mysound(male_short,Fs);
% display_(male_short,Fs,1,0);

figure, subplot(2,1,1)
spectrogram(S, 200, 50, 200, Fs, 'yaxis'); colorbar('off');
title(['male short Frequency vs time']);
subplot(2,1,2);
display_(S, Fs, scaling_f, mute);
%%

%% DFT
% x is a vector of speech, N is the frame length (must be even),
% S is the first sample of the frame to be analyzed


S=1;
x = male_short;
N = 100;
xf = x(S:S+N-1).*hanning(N);
X = fft(xf);
figure(1); clf;
plot(10*log10(abs(X(1:N/2+1)).^2));