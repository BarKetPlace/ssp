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

    % loop over all windows but the last one which is always smaller than
    % the rest therefore cannot be put inside the matrix S
for n= 1 : numwin - 1
        % select part of signal
    xf = x(n1:n2) ;
    
        % compute DFT
    X = fft(xf.*hanning(winlen));
    S(:, n) = 10*log10(abs(X(1:winlen / 2 + 1).^2)) ;
    
        % show DFT + time signal of window
    figure(1), clf
    subplot(211)
    plot(xf)
    title('Time plot of window')
    axis([1 winlen min(x) max(x)])
    
    subplot(212)
    plot(10*log10(abs(X.^2)));
    axis([1 winlen -25 125])
    title('DFT of window')    
    
    pause(0.01)
    
        % update boundaries of window
    n1 = n1 + uplen;
    n2 = n2 + uplen;
end

    % plot spectrogramme
figure,
colormap(gray)
imagesc(flipud(-S));




