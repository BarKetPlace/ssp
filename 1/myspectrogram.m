function S = myspectrogram(x, Fs, winlen, uplen)
%MYSPECTROGAM Computes and plot the spectrogram of signal x
%   x       : Signal to analyse
%   Fs      : Sampling frequency of signal
%   winlen  : Size of window analysis
%   uplen   : Number of updated samples

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
    
        % compute DFT and store it in matrix
    X = fft(xf.*hanning(winlen), N);
    S(:, n) = 10*log10(abs(X(winlen/2:end).^2))' ;
    
        % update boundaries of window
    n1 = n1 + uplen;
    n2 = n2 + uplen;
end

    % plot spectrogramme
figure,
colormap(bone)
imagesc([0 1000*N/Fs], [Fs/2 0], flipud(S));
c = colorbar;
c.Label.String = 'Values (dB)';
title('Spectrogram')
ylabel('Frequency (Hz)')
xlabel('Time (ms)')


end

