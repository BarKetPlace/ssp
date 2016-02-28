clear all
close all
clc

load assignment2.mat

in = [-6:.01:6]'; %Signal must be a columns
m=0; %Offset of the quantizer recosntruction levels =0 ->midrise quantizer =delta/2 ->midtread
xmax = 4;
n_bits = 2;
en_plots = 1; %Enables plot in functions

[ idx ] = sq_enc(in, n_bits, xmax, m, en_plots);

outq = sq_dec(idx, n_bits, xmax, m);


figure, plot(in,outq);
title(['Uniform Scalar Quantizer. ' num2str(n_bits) ' bits, m= ' num2str(m) ', xmax= ' num2str(xmax) ]);
xlabel('Input');ylabel('Output');
