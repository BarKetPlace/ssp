function [E, ZC, V, A, P] = analysis(x, alen, ulen, M, Fs)

    % get length
N = length(x) ;            % samples
ms = (N / Fs) * 1000 ;       % milliseconde

    % get size of last window
lastwin = rem( N - alen, ulen) ;

    % get number of windows
naf = (N - alen - lastwin) / ulen + 2 ;
naf = naf - 1 ;   % the last window is not taken into account because not the same size
    
    % Initialization of variables
E = zeros(naf, 1);      % energy
ZC = zeros(naf, 1);     % zerocrossings number
V = zeros(naf, 1);      % voiced/unvoiced
A = zeros(naf, M+1);    % LP analysis
P = zeros(naf,1);       % pitch
voicedthreshold = .3 ;  % threshold for voiced/unvoiced detection

    % initialize boundaries of window
n1 = 1 ;
n2 = alen ;

for n = 1 : naf
        % select part of signal
    xf = x(n1:n2) ;

        %%% 1- get energy
    E(n) = 1/alen*sum(xf.^2);
    
        %%% 2- get voiced/unvoiced detection
    for i = 1 : alen -1
        if(xf(i)*xf(i+1) < 0 )
                % increment ZC each time we have a zerocrossing
            ZC(n) = ZC(n) + 1 ; 
        end
    end
    ZC(n) = ZC(n)/alen;                 % Normalization
    V(n) = ZC(n) <= voicedthreshold ;   % Decision making
    
        %%% 3- get LP analysis
    xf_win = xf .* window(@hanning, alen) ;     % windowing of xf
    c = xcorr(xf_win, xf_win ,M) ;              % get autocorrelation
    [a, ~]= levinson( c(M+1:2*M+1) ) ;          % get polynomial coefficients
    A(n,:) = a ;                                % store them in A
    
        %%% 4- get pitch with ACF
            % get half of correlation since correlation is symmetric
    c = xcorr(xf, xf) ;
    c = c(round(alen-1/2):end) ;
            % get all the peaks
    peaks = logical([0; (c(2:end-1) > c(1:end-2)) & (c(2:end-1) > c(3:end)); 0]) ;
            % find the highest peak
    ind = find(c == max(c(peaks)) ) ;
            % compute the pitch
    P(n) = Fs / ind ;       % Hz
    
            % plot for accuracy
%     plot(c)
%     line([ind ind], [min(c) max(c)], 'Color', 'red') ;
%     pause(0.001)
    
        % update boundaries of window
    n1 = n1 + ulen;
    n2 = n2 + ulen;
end

    % plot results
x1 = linspace(0, ms, naf) ;
figure(1);clf;

        % Plot the input waveform
subplot(3,2,1)
plot(linspace(0, ms, N), x)
axis([1 ms min(x) max(x)]);
title('Time plot of signal');
xlabel('time (ms)')

        % Plot the standard deviation
subplot(3,2,2)
plot(x1, sqrt(E)) 
axis([1 ms min(sqrt(E)) max(sqrt(E))]);
title('Standard deviation');
xlabel('time (ms)')

        % Plot voiced/unvoiced decision
subplot(3,2,3)
plot(x1, V) 
axis([1 ms 0 1.5]);
title('Voiced (1) unvoiced(0)');
xlabel('time (ms)')

        % Plot the normalized number of zero-crossings
subplot(3,2,4)
plot(x1, ZC); hold on;
plot(length(ZC)*[0 1], voicedthreshold*[1 1], '-r');
axis([1 ms min(ZC) max(ZC)]);
title('Normalized number of ZC');
xlabel('time (ms)')

        % Plot the fundamental frequency in Hz
subplot(3,2,5)
plot(x1, P) 
axis([1 ms min(P) max(P)])
ylabel('Frequency (Hz)');
title('Fundamental frequency');
xlabel('time (ms)')

        % Illustrate the vocal tract envelope in a spectrogram style!
subplot(3,2,6)
S = zeros(512, naf);
for n=1:naf
S(:,n) = 20*log10(abs(freqz(1,A(n,:),512)));
end
S = flipud(S);
colormap(gray);
imagesc(S); 
title('Vocal tract envelope');

end