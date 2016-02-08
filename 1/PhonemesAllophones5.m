close all
clear all
clc

load assignment1.mat

S = male_short;
Fs = 8000;
scaling_f = 1;
mute = 0;

part=1:length(S);
% part=700:1800;
figure, subplot(2,1,1)
spectrogram(S(part), 200, 50, 200, Fs, 'yaxis'); colorbar('off');
title(['male short Frequency vs time']);
subplot(2,1,2);
display_(S(part), Fs, scaling_f, mute); hold on;
min_ = min(S); max_ = max(S);

figure, plot(S)
figure,
plot(S);hold on;
%[n]
plot(700*[1 1], [min_ max_],'-r'); hold on;
%[ow]
plot(1800*[1 1], [min_ max_],'-r'); hold on; 
%[he's]
plot(2800*[1 1], [min_ max_],'-r'); hold on;
%[oua]
plot(3300*[1 1], [min_ max_],'-r'); hold on; 
%[n]
plot(3700*[1 1], [min_ max_],'-r'); hold on; 
%[o]
plot(3900*[1 1], [min_ max_],'-r'); hold on; 
%[f]
plot(4600*[1 1], [min_ max_],'-r'); hold on; 
%[the]
plot(5200*[1 1], [min_ max_],'-r'); hold on; 
%[m]
plot(5600*[1 1], [min_ max_],'-r'); hold on; 
%[o]
plot(6500*[1 1], [min_ max_],'-r'); hold on;
%[s]
plot(7000*[1 1], [min_ max_],'-r'); hold on; 
%%
% x is a vector of speech, N is the frame length (must be even),
% S is the first sample of the frame to be analyzed

%Voiced : [v][g][b][z][d][>th<is]
%Sound n from the word "Now"
type_{1} = 'Voiced';
S_(1)=122;
N_(1) =450;
%Unvoiced : [f] [k] [p][s][t][tch][th][ss]
% Sound s from the word "moSt"
type_{2} = 'Unvoiced';
S_(2)=5900;
N_(2) = 1000;
x = male_short;