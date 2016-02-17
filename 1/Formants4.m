%% load stuff
clear all
close all
clc

    % load vectors
load assignment1.mat
load fig1formants.mat

    % set sampling frequency
Fs = 8000 ;     % Hz

%% Perform LP analysis
close all
clc

    %%%% PARAMETERS TO SET %%%%
signal = male_short ;                   % choose the signal to analyse
N = 200 ;                              % window length
M = 10  ;                               % prediction order

    %%%% OTHER PARAMETERS %%%%
lengthsignal = length(male_short) ;     % length of signal
lastwind =  rem(lengthsignal, N) ;      % length of last window
Nwin = (lengthsignal - lastwind) / N ;  % number of windows
indices = 1:N:lengthsignal ;            % indices of first sample for each frame
vowel = [] ;                            % cell of vowels
j = 1 ;                                 % indice of windows

    % for each window, we want to find the formant
for i = indices
    
        % if we're at the last window, which is smaller, we have to change
        % the length to lastwind
    S = N ;
    if i == indices(end)
        S = lastwind ;
    end
    
        % get windowed portion of signal
    xf = signal(i:i+S-1).*window(@hanning, S) ;
    
        % get autocorrelation
    c = xcorr(xf, xf, M) ;
        
        % get coefficients ai for LP analysis
    [a, e]= levinson(c(M+1:2*M+1)); % a is a vector always starting with a 1.
    a = a(:); % Make a a column vector

    h = freqz(1, a, S, 'whole');
    
        % plot result
    figure
    plot(linspace(0, Fs/2, length(1:S/2+1)), 10*log10(e*abs(h(1:S/2+1)).^2), 'LineWidth',1.5) ;
    title('Use the cursor to pick the 3 formants of this window')
     
        % for each plot, the user is asked to find formants and answer with
        % vowel
    fprintf('Window %d / %d\n', j, Nwin + 1) ;
    [x, ~]  = ginput(3) ;
    
        % select which vowel is the closest
    mindist = Inf ;
    for k=1:length(fig1formants)
   
            % get formant from structure
        formant = [fig1formants(k).F1; fig1formants(k).F2; fig1formants(k).F3] ;

            % get closesness to that vowel
        a = sum((x-formant).^2) ;

        if a < mindist
            mindist = a ;
            indformant = k ;
        end    
    end

        % get vowel corresponding to the closest formant
    vowel = [vowel; fig1formants(indformant).vowels] ;
    
        % close window and increment j
    close all
    j = j + 1;    
end


%% time plot signal with vowels
close all
clc

Ymin = -10000 ;
Ymax = 15000 ;
Ytext = 14000 ;

figure
plot(1:lengthsignal, signal)
title(['Time plot of signal and vowel analysis with window length of ', int2str(N)])
xlabel('Time (sample)')
axis([1 lengthsignal Ymin Ymax])

    % draw vertical lines that bound the windowed samples
for i = indices
    line([i, i], [Ymin, Ymax], 'LineWidth', 1, 'LineStyle', ':', 'Color', 'black') ;   
end

    % add vowels
X = indices + N/2 ;
X(end) = (indices(end) + lengthsignal)/2 ;
text(X, Ytext*ones(Nwin+1, 1), vowel) ;

