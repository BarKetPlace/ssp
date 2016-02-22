%% Load stuff
clear all
close all
clc

load assignment1.mat
[female,Fss] = audioread('female44.wav');     % female sound
male = audioread('male44.wav');         % male sound
Fs = 8000 ;                             % sample frequency of signal
soundsc(female,Fss);
%% 2- Bandwith of Speech
close all
clc
fprintf('Cut off frequency     \n');
for i=1200:300:2000
        % set cut off frequency
    Fc = i / Fs ;      % normalized cut off frequency

        % filter the signal
    filtered_ = lowpass(female, Fc) ;
    soundsc(filtered_,Fss);
    fprintf('\b\b\b\b\b%04d\n',i);
    pause(length(filtered_)/Fss);
        % write file
%     audiowrite(['audio/filtered_male_', int2str(i), '.wav'], filtered_male, Fss) ;
end

fprintf('The signal is filtered with a cutoff frequency from 100 to 500 Hz with step of 50 Hz\n') ;
fprintf('We start to distinguish words from 250 Hz but the full message is clear at 400 Hz\n')
% 1) 3000Hz we understand the message perfectly
% 2) 400Hz the message is not understandable
% 3) The signal is perfectly understandable with a cut off frequency of
% 3000Hz. The rule is that Fs>= 2*Fm where Fm is the highest frequency of
% the signal