close all
clear all
clc

load assignment2.mat

Fs = 8000;
scaling_f = 1;
mute = 0;
	% set parameters
x = male8 ;
flen=30 ;               %ms
alen = flen/1000*Fs ;   % size of analysis window
ulen = 120 ;             % length of update
M = 120;                 % order for LP analysis

en_plots=1;
% soundsc(x,Fs);
[E, V, A,P]=analysis(x, alen, ulen, M);%, Fs);%, en_plots);
s = synthesis(E, V, A, P, ulen);
% s = s/max(s)*max(x);

figure, 
plot(x); hold on;
plot(s); title(['ulen = ' num2str(ulen)]);
legend('Initial', 'Vocoded');
soundsc(s,Fs);