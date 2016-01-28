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

for i=100:50:500
        % set cut off frequency
    Fc = i / Fs ;      % normalized cut off frequency

        % filter the signal
    filtered_male = lowpass(male, Fc) ;
    
        % write file
    audiowrite(['audio/filtered_male_', int2str(i), '.wav'], filtered_male, Fss) ;
end

fprintf('The signal is filtered with a cutoff frequency from 100 to 500 Hz with step of 50 Hz\n') ;
fprintf('We start to distinguish words from 250 Hz but the full message is clear at 400 Hz\n')

%% 
