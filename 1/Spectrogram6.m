%% CREATE A SPECTOGRAM FUNCTION
clear all
close all
clc

load assignment1.mat ;

%% ANSWER TO QUESTIONS 1) AND 2)
close all
clc

    % get input variables
x = male_short ;            % signal
winlen = 256 ;              % window length
uplen = 32 ;               % overlap factor [0,1] => percentage of winlen

    % get length
N = length(x) ;

    % get size of last window
lastwin = rem( N - winlen, uplen) ;

    % get number of windows
numwin = (N - winlen - lastwin) / uplen + 2 ;

    % create spectogram matrix S
S = zeros(winlen / 2 + 1, numwin) ;

    % initialize boundaries of window
n1 = 1;
n2 = winlen;

    % set parameters for plot
timewin = winlen/(1000*Fs) ;            % ms
Xtim = linspace(0, timewin, winlen) ;   % abscissa for time plot
Xfreq = linspace(-Fs/2, Fs/2, winlen) ; % abscissa for freq plot

    % loop over all windows but the last one which is always smaller than
    % the rest therefore cannot be put inside the matrix S
for n= 1 : numwin - 1
        % select part of signal
    xf = x(n1:n2) ;
    
        % compute DFT
    X = fft(xf.*hanning(winlen));
    S(:, n) = 10*log10(abs(X(winlen/2:end).^2))' ;
    
        % show DFT + time signal of window
    figure(1), clf
    subplot(211)
    plot(Xtim, xf)
    title('Time plot of window')
    xlabel('Time (ms)')
    axis([0 timewin min(x) max(x)])
    
    subplot(212)
    plot(Xfreq, 10*log10(abs(X.^2)));
    axis([-Fs/2 Fs/2 -25 125])
    title('DFT of window')
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dB)')
    
    pause(0.01)
    
        % update boundaries of window
    n1 = n1 + uplen;
    n2 = n2 + uplen;
end

%% spectogramme
close all

    % plot spectrogramme
figure,
colormap(gray)
imagesc([0 N / (1000*Fs)], [0 Fs/2], flipud(S));
c = colorbar;
c.Label.String = 'Values (dB)';
title('Spectrogramme')
ylabel('Frequency (Hz)')
xlabel('Time (ms)')


%% test of function
close all
clc

bits = 12 ;

myspectrogram(male_short, Fs, 2^bits, 2^(bits-3)) ;




