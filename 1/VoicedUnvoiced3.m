%%3- Voiced and Unvoiced Speech Sounds
clear all
close all
clc

load assignment1.mat

S = male_short;
Fs = 8000;
scaling_f = Fs;
mute = 0;

figure, subplot(2,1,1)
spectrogram(S, 200, 50, 200, Fs, 'yaxis'); colorbar('off');
title(['male short Frequency vs time']);
subplot(2,1,2);
display_(S, Fs, scaling_f, mute);
%%

%% DFT
close all

% x is a vector of speech, N is the frame length (must be even),
% S is the first sample of the frame to be analyzed

%Voiced : [v][g][b][z][d][>th<is]
%Sound n from the word "Now"
type_{1} = 'Voiced';
S_(1)=122;
N_(1) =450;
%Unvoiced : [f] [k] [p][s][t][tch][th][ss]
% Sound s from the word "moSt"
type_{2} = 'Unvoiced';
S_(2)=5900;
N_(2) = 1000;
x = male_short;
for itype = 1:2
    S = S_(itype);
    N = N_(itype);
    
    xf = x(S:S+N-1).*hanning(N);
    
    X = fft(xf);
    
%     axis_x = [0:1/scaling_f:(N-1)/scaling_f];
    axis_x = [0:2/N:1].*Fs/2;
%     axis_x = axis_x(1:end/2+1);
    % 3) PLOT
    figure,
    plot(axis_x,10*log10(abs(X(1:N/2+1)).^2),'LineWidth',1); hold on;
    % plot(10*log10(abs(Xv_rect(1:Nv/2+1)).^2),'LineWidth',1.5);
    
    xlabel('Frequency(Hz)');
    ylabel('Amplitude (dB)');
    title([type_(itype) 'frame length : ' num2str(length(xf)) 'Hanning window']);
    
    
    % 4) k=3 <> 2*pi*k/Nv*8000 = 335
    k=3; 2*pi*9/N*8000;
    % 5) First pic k=9 <> pitch is 1000Hz
    % 6) Herethe frame length is 450. If the frame length is too short, the signal
    % will not have time to end one period and the result will be irrelevant.
    % If the frame length is too short, the ergodic approximation will not hold
    % 7) More power with the hanning window
    
    M_ = [1, 10, 100];
    for indM = 1:length(M_)
        M = M_(indM);
        legendM{indM} = ['Enveloppe, Prediction order :: ' num2str(M_(indM))];
        % M is the prediction order
        c = xcorr(xf, xf, M);
        [a, e]= levinson(c(M+1:2*M+1)); % a is a vector always starting
        %with a 1.
        a = a(:); % Make a a column vector
        h = zeros(N,1);
        for k=0:N-1
            h(k+1) = 1 / (a'*exp(-i*2*pi*k/N*(0:M)') );
        end
        % h = freqz(1, a, N, ’whole’); % Equivalent to the above for-loop!
        hold on
        plot(axis_x,10*log10(e*abs(h(1:N/2+1)).^2), 'LineWidth',1.5);%, 'r');
    end
    
    legend('Spectrum',legendM{1:3});
end