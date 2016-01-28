%% Load stuff
clear all
close all
clc

load assignment1.mat
female = audioread('female44.wav');     % female sound
male = audioread('male44.wav');         % male sound
Fs = 8000 ;                             % sample frequency of signal
Fss = 44100 ;                           % sample frequency of speech

%% 2- Bandwith of Speech
close all
clc

for Fc=100:50:500
        % set cut off frequency
%     Fc = 100 ;          % Hz
    Fc = Fc / Fs ;      % normalized cut off frequency

        % filter the signal
    filtered_female = lowpass(female, Fc) ;

        % listen to the signal
    soundsc(filtered_female, Fss) ;

   pause(3)
end


ssadasd
