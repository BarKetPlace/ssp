clear all
close all
clc

load assignment2.mat

x = male8;
Fs = 8000;
N = 30*Fs/1000; %Analysis window length
U = N; %Update window length
M = 10; % Orderfor linear prediction

en_plots = 1;

%Compute linear prediction error
n_frames = floor((length(x)-N+U)/U); % Number of frame
ms = (length(x) / Fs) * 1000 ;       % milliseconde
x1 = linspace(0, ms, n_frames) ;
A = zeros(n_frames, M+1);
Err = zeros(size(x));
s = 1; % Index for the start of the analysis frame
e = N;% Index for the end of the analysis frame
% Make x a column vector
x = x(:);
for n=1:n_frames
    xf = x(s:e);
    
    % Linear prediction analysis:
    a = lpc(xf,M); %The a coeff are directly the openloop system coefficients
    a = real(a);
    A(n,:) = a;
    % end lp analysis
    
    Err(s:e,1) = filter(A(n,:),1,xf); %Filter with 1-A(Z)
%     Err(s:e) = xf - y(s:e);
    s = s + N;
    e = e + N;
end
%% Plot figure, and listen
sigcorr = xcorr(x);
errcorr = xcorr(Err);

figure, %Observe the decorrelation
plot(sigcorr/max(sigcorr)); hold on;
plot(errcorr/max(errcorr), 'LineWidth', 2);
ylabel('Auto-Correlation values');
xlabel('Delay (in samples)'); 
title('Auto-correlation function of the speech signal and error signal');
legend('Speech signal','Error signal');

soundsc(Err,Fs);
%% Encoding of A and Err


%% Decoding of A and Err


%% Receiver
s = 1; % Index for the start of the synthesis frame
e = U;% Index for the end of the synthesis frame
% Make Err a column vector
Err = Err(:);

for n=1:n_frames
    Errf = Err(s:e);
    
    xhat(s:e) = filter(1,A(n,:), Errf);
    %y(s:e,1) = filter([1 -a(2:end)],1,xf);
%     Err(s:e) = xf - y(s:e);
    s = s + U;
    e = e + U;
end
figure, plot(xhat); hold on;
    plot(x);
    legend('xhat','x');

  soundsc(xhat,Fs);
  