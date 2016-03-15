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
M    =  10;                            % order for LP analysis
ms   =  (length(x) / Fs) * 1000 ;       % milliseconde
D   =   [] ;                            % distortion 

    % Analysis
[E, V, A, P]=analysis(x, Fs, alen, ulen, M, 1) ;

    %% 3.1 - Quantize gain as it is
close all
clc

    % create figures
% fig1 = figure('Position',[1 1 scrsz(3)/2 scrsz(4)]) ;
fig1 = figure
fig2 = figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]) ;


% for nbits = 8:-1:2    
        % display bits
    fprintf(['Bits = ', num2str(nbits), '\n'])
    
        % plot histogram
    figure(fig1) ;
    a = histogram(E) ;
        % get m and xmax
    m = mean(a.BinEdges) ;
    xmax = max(a.BinEdges) - m ;
        % plot vertical lines
    hold on
    plot([m-xmax, m-xmax], [0 max(a.Values)], 'red', [m, m], [0 max(a.Values)], 'red', [m+xmax, m+xmax], [0 max(a.Values)], 'red') ;
    text(m-xmax+5, max(a.Values)-5, 'm-xmax'), text(m+5, max(a.Values)-5, 'm'), text(m+xmax+5, max(a.Values)-5, 'm+xmax')
    title('Histogram of Gain')
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
        
% end

    % result
fprintf('Can''t hear a difference after 4 bits\n')


    %% 3.1 - Quantize log of gain
close all
clc

    % create figures
% fig1 = figure('Position',[1 1 scrsz(3)/2 scrsz(4)]) ;
% fig2 = figure('Position',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)]) ;
fig1 = figure,
fig2 = figure,

for nbits_E = 8:-1:3  
        % display bits
    fprintf(['Bits = ', num2str(nbits_E), '\n'])
    
        % plot histogram
    figure(fig1) ;
    a = histogram(log(E)) ;
        % get m and xmax
    m = mean(a.BinEdges) ;
    xmax = max(a.BinEdges)-m ;
    
        % plot vertical lines
    hold on
    plot([m-xmax, m-xmax], [0 max(a.Values)], 'red', [m, m], [0 max(a.Values)], 'red', [m+xmax, m+xmax], [0 max(a.Values)], 'red') ;
    text(m-xmax, max(a.Values), 'm-xmax'), text(m, max(a.Values), 'm'), text(m+xmax, max(a.Values), 'm+xmax')
    title('Histogram of the log of the gain')
    hold off

        % quantize
    idx = sq_enc(log(E), nbits_E, xmax, m, 0) ;
    quantlogE = sq_dec(idx, nbits_E, xmax, m) ;
    
        % normal quantization
    a = histogram(E) ;
        % get m and xmax
    m = mean(a.BinEdges) ;
    xmax = max(a.BinEdges) - m ;
            % quantize
    idx = sq_enc(E, nbits_E, xmax, m, 0) ;
    quantE = sq_dec(idx, nbits_E, xmax, m) ;

        % plot the result and compare to quantization values and original
        % signal
    figure(fig2)
    plot(E); hold on; plot(exp(quantlogE)) ; plot(quantE) ; hold off ;
    legend('Gain', 'Log quantization of gain', 'Direct quantization of gain')
    title(['Quantizations of gain with bits = ', num2str(nbits_E)])
    
        % synthesis with parameters
%     x_synth = synthesis(exp(quantE), V, A, P, ulen) ;
%     
%         % play
%     soundsc(x_synth, Fs) ;
%     pause(ms/1000)
        
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
    
%%%%%% CODE FUNCTION

    % create codeA vector
codeA = zeros(length(A), 2) ;

    % for each polynomial vectors
for i = 1:length(A)
   
        % get lsf
    lsf = poly2lsf(A(i,:))' ;
    
        % get distance between lsf and codebook
    dist = (1/M) * sum((ones(2^10, 1)*lsf - lsfCB1).^2, 2) ;
    
        % get index of smallest distance
    indx1 = find( dist == min(dist)) ;
    
        % get residual and code it same way as before but with lsfCB2
    res = lsf - lsfCB1(indx1, : ) ;
    dist = (1/M) * sum((ones(2^10, 1)*res - lsfCB2).^2, 2) ;
    indx2 = find( dist == min(dist)) ;
    
        % update codeA vector with indexes from both codebooks
    codeA(i, :) = [indx1, indx2] ;        
end

%%%%% DECODE FUNCTION

    % create vector Aq
Aq = zeros(size(A)) ;

for i = 1:length(codeA)
   
        % get indexes
    indxs = codeA(i, :) ;
    
        % get lsf et residual
    lsf = lsfCB1(indxs(1), :) ;
    res = lsfCB2(indxs(2), :) ;
    lsf = sort(lsf + res) ;         % use sort to make sure it represents a minimum phase whitening filter

        % call lsf2poly and store the result
    Aq(i, :) = lsf2poly(lsf) ;
end

    %% 3.4 - Optimizing the bit allocation
close all
clc

    % set quantization parameters
nbits_E = 8 ;                       % 2 bits for gain
nbits_P = 8 ;                       % 4 bits for pitch
SNR = zeros(nbits_E, nbits_P) ;
for nbits_E=1:8
    for nbits_P=1:8
nbits_V = 1 ;                       % 1 bit for Voiced/Unvoiced
nbits_LP = 2 * M ;                  % fixed number of bits for LP parameters quantization
                                    % depends on number of window and LP order

    % get rate
nbits_total         = length(A)*(nbits_E + nbits_P + nbits_V + nbits_LP) ;
ratebitpersample    = nbits_total / length(x) ;                     % bit   / sample 
ratebitperseconde   = nbits_total / (ms) ;                          % kbits / sec

    % GAIN QUANTIZATION
a = histogram(log(E), 2^nbits_E) ;
m = mean(a.BinEdges) ;
xmax = max(a.BinEdges)-m ;
idx = sq_enc(log(E), nbits_E, xmax, m, 0) ;
Eq = sq_dec(idx, nbits_E, xmax, m) ;
Eq = exp(Eq) ;

    % VOICED/UNVOICED QUANTIZATION
Vq = V ;

    % PITCH QUANTIZATION
a = histogram(P, 2^nbits_P) ;
m = mean(a.BinEdges) ; 
xmax = max(a.BinEdges)-m ; close all
idx = sq_enc(P, nbits_P, xmax, m, 0) ;
Pq = sq_dec(idx, nbits_P, xmax, m) ;
% Pq = exp(Pq) ;

    % LP PARAMETERS QUANTIZATION
codeA = encodefilter(A, lsfCB1, lsfCB2) ;
Aq = decodefilter(codeA, lsfCB1, lsfCB2) ;

    % SYNTHESIS
outputx = synthesis(Eq, Vq, A, Pq , ulen) ;

    % PLOT
delay = length(x) - length(outputx) ;
figure,
plot(x) ; hold on;
plot(outputx) ;
% soundsc(outputx, Fs) ; pause(ms/1000) ;
% soundsc(x, Fs) ;

SNR(nbits_E,nbits_P) = 10*log10(var(x)/var(outputx)) ;
    end
end

%%

    % compute SNR
% SNR = (1/length(x)) * sum((x(delay+1:end)-outputx).^2) ;
SNR = zeros(delay+1,1) ;
% figure
for i= 0:delay
%     plot(x(1+i:end-delay+i)) ; hold on;
%     plot(outputx) ;
%     axis([1000 2300 -8000 10000])
%     hold off;
%     pause(0.1)

    SNR(i+1) = 10*log10((1/length(x)) * sum((x(1+i:end-delay+i)-outputx).^2)) ;
      
end
    