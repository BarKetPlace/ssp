clear all
close all
clc


load assignment1.mat

S = male_short;
Fs = 8000;
scaling_f = 1;
mute = 0;

x = S;
%Frame length 
flen=30;%ms
alen = 256;% flen/1000*Fs;
voicedthreshold = .4; %Thresholdvoiced/unvoiced
naf = ceil(length(x)/alen);

ulen =32;%20;
M = 1;
% mysound(x,Fs);
[E, ZC, V, A,P]=analysis(x,alen,ulen,voicedthreshold,M);

s = synthesis1(E, ZC, V, A, P, ulen);
figure, plot(s); title(['ulen = ' num2str(ulen)]);
mysound(s,Fs);


