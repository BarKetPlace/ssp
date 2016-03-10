    %% load stuff
close all
clear all
clc

load assignment2.mat                % load vectors
Fs = 8000 ;                         % sampling frequency
scrsz = get(groot,'ScreenSize') ;   % get screen size

    %% Perform analysis to get Gain, Pitch, Voiced/Unvoiced and LP
close all
clc

	% set parameters
x    =  male8 ;
flen =  30 ;                            % ms
alen =  flen/1000*Fs ;                  % size of analysis window
ulen =  120 ;                           % length of update
M    =  120;                            % order for LP analysis
ms   =  (length(x) / Fs) * 1000 ;       % milliseconde
D   =   [] ;                            % distortion 

    % Analysis
[E, V, A, P]=analysis(x, Fs, alen, ulen, M, 1) ;

    %% 3.1 - Quantize gain as it is
close all
clc

    % create figures
fig1 = figure('Position',[1 1 scrsz(3)/2 scrsz(4)]) ;
fig2 = figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]) ;

for nbits = 8:-1:2    
        % display bits
    fprintf(['Bits = ', num2str(nbits), '\n'])
    
        % plot histogram
    figure(fig1) ;
    a = histogram(E, 2^nbits) ;
        % get m and xmax
    m = mean(a.BinEdges) ;
    xmax = max(a.BinEdges) - m ;
        % plot vertical lines
    hold on
    plot([m-xmax, m-xmax], [0 max(a.Values)], 'red', [m, m], [0 max(a.Values)], 'red', [m+xmax, m+xmax], [0 max(a.Values)], 'red') ;
    text(m-xmax+5, max(a.Values)-5, 'm-xmax'), text(m+5, max(a.Values)-5, 'm'), text(m+xmax+5, max(a.Values)-5, 'm+xmax')
    hold off

        % quantize
    idx = sq_enc(E, nbits, xmax, m, 0) ;
    quantE = sq_dec(idx, nbits, xmax, m) ;

        % plot the result and compare to quantization values and original
        % signal
    figure(fig2)
    plot(E); hold on; plot(quantE) ; hold off ;
    legend('Gain', 'Quantized gain')
    title(['Quantization of gain with bits = ', num2str(nbits)])

        % synthesis with parameters
    x_synth = synthesis(quantE, V, A, P, ulen) ;
    
        % play
    soundsc(x_synth, Fs) ;
    pause(ms/1000)
        
end

    % result
fprintf('Can''t hear a difference after 4 bits\n')


    %% 3.1 - Quantize log of gain
close all
clc

    % create figures
fig1 = figure('Position',[1 1 scrsz(3)/2 scrsz(4)]) ;
fig2 = figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]) ;

for nbits_E = 8:-1:2    
        % display bits
    fprintf(['Bits = ', num2str(nbits_E), '\n'])
    
        % plot histogram
    figure(fig1) ;
    a = histogram(log(E), 2^nbits_E) ;
        % get m and xmax
    m = mean(a.BinEdges) ;
    xmax = max(a.BinEdges)-m ;
        % plot vertical lines
    hold on
    plot([m-xmax, m-xmax], [0 max(a.Values)], 'red', [m, m], [0 max(a.Values)], 'red', [m+xmax, m+xmax], [0 max(a.Values)], 'red') ;
    text(m-xmax, max(a.Values), 'm-xmax'), text(m, max(a.Values), 'm'), text(m+xmax, max(a.Values), 'm+xmax')
    hold off

        % quantize
    idx = sq_enc(log(E), nbits_E, xmax, m, 0) ;
    quantE = sq_dec(idx, nbits_E, xmax, m) ;

        % plot the result and compare to quantization values and original
        % signal
    figure(fig2)
    plot(E); hold on; plot(exp(quantE)) ; hold off ;
    legend('Gain', 'Quantized gain')
    title(['Quantization of log of gain with bits = ', num2str(nbits_E)])
    
        % synthesis with parameters
    x_synth = synthesis(exp(quantE), V, A, P, ulen) ;
    
        % play
    soundsc(x_synth, Fs) ;
    pause(ms/1000)
        
end

    % result
fprintf('Can''t hear a difference after 4 bits\n') ;


%%%%% Commentary
% 
% better with log
% plan : using log quantization with 3 bits

    %% 3.2 - Quantize Pitch & Voiced/Unvoiced
close all
clc

%%%% Voiced/Unvoiced is already a 1 bit signal since it represents a
%%%% binary decision, we can't quantize it more without loosing crucial
%%%% information

nbits_V = 1 ;

%%%% in order to quantize the Pitch correctly we are going to consider the
%%%% logarithmic human perception of sounds

   % This vector contains the frequencies of the scales 0 to 4
frequencies = (-57:14)/12  ;
x = 2*ones(1,length(frequencies)) ;
frequencies = (x.^frequencies)*440 ;
          

    % create figures
fig1 = figure('Position',[1 1 scrsz(3)/2 scrsz(4)]) ;
fig2 = figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]) ;

for nbits_P = 8:-1:1    
        % display bits
    fprintf(['Bits = ', num2str(nbits_P), '\n'])
    
        % filter P
    
        % plot histogram
    figure(fig1) ;
    a = histogram(log(P), 2^nbits_P) ;
        % get m and xmax
    m = mean(a.BinEdges) ;
    xmax = max(a.BinEdges)-m ;
        % plot vertical lines
    hold on
    plot([m-xmax, m-xmax], [0 max(a.Values)], 'red', [m, m], [0 max(a.Values)], 'red', [m+xmax, m+xmax], [0 max(a.Values)], 'red') ;
    text(m-xmax, max(a.Values), 'm-xmax'), text(m, max(a.Values), 'm'), text(m+xmax, max(a.Values), 'm+xmax')
    hold off

        % quantize
    idx = sq_enc(log(P), nbits_P, xmax, m, 0) ;
    quantP = sq_dec(idx, nbits_P, xmax, m) ;

        % plot the result and compare to quantization values and original
        % signal
    figure(fig2)
    plot(log(P)); hold on; plot(quantP) ; hold off ;
    legend('Gain', 'Quantized gain')
    title(['Quantization of log of pitch with bits = ', num2str(nbits_P)])
    
        % synthesis with parameters
    x_synth = synthesis(E, V, A, exp(quantP) , ulen) ;
    
        % play
    soundsc(x_synth, Fs) ;
    pause(ms/1000)
        
end

% we notice a clear change in pitch under 4 bits quantization, it seems
% that under 16 levels of quantization we deteriorate the signal too much

    %% 3.3 - Quantizing the LP parameters
    
    