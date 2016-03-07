close all
clear all
clc

    % load stuff
load assignment2.mat
Fs = 8000 ;     % sampling frequency

    %% quantize gain
close all
clc

	% set parameters
x    =  male8 ;
flen =  30 ;                            % ms
alen =  flen/1000*Fs ;                  % size of analysis window
ulen =  120 ;                           % length of update
M    =  120;                            % order for LP analysis
ms   =  (length(x) / Fs) * 1000 ;       % milliseconde

    % Analysis
[E, V, A, P]=analysis(x, Fs, alen, ulen, M, 0) ;

    % plot gain
plot(linspace(0, ms, length(E)), E) ;
axis([1 ms min(E) max(E)]);
title('Gain');
xlabel('time (ms)')

    % plot histogram of gain
nbits = 5 ;
figure
a = histogram(E, 2^nbits);

