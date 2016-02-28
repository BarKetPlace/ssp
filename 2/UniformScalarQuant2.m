clear all
close all
clc

load assignment2.mat

in = male8; %Signal
m=0; %Mean of the possible values 
xmax = max(male8);
n_bits = 3;
en_plots = 0; %Enables plotsin functions /!\ it is very expensive /!\

[ idx ] = sq_enc(in, n_bits, xmax, m, en_plots);

outq = sq_dec(idx, n_bits, xmax, m);


figure, plot(in); hold on; plot(outq);
title('Quantization result');
legend('USQ input', 'USQ output');
