clear all
close all
clc


load assignment1.mat
Fs = 8000;
scaling_f = 1;
mute = 0;
	% set parameters
x = male_short ;
flen=30 ;               %ms
alen = flen/1000*Fs ;   % size of analysis window
ulen = 70 ;             % length of update
M = 50;                 % order for LP analysis

en_plots=1;
mysound(x,Fs);
[E, ZC, V, A,P]=analysis(x, alen, ulen, M, Fs, en_plots);
s = synthesis2(E, ZC, V, A, P, ulen);
% s = s/max(s)*max(x);

figure, 
plot(x); hold on;
plot(s) title(['ulen = ' num2str(ulen)]);
legend('Initial', 'Vocoded');
mysound(s,Fs);


