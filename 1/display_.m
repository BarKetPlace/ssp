%function [] = display_(S,Fs,scaling_f,mute)
%display_ 
% subplot(2,1) ::   upper plot is the spectrogram(frequency vs time)
%                   lower plot is the temporal signal (with a player)
%
% S is a vector containing the samples to play
% Fs is the sampling frequency
% scaling_f = Fs to display the x axis in sec, 
%           = 1  to display x axis in samples
%%----------------------------------------------------
%Code Authors:
% Antoine Honoré
% Audrey Brouard
%----------------------------------------------------

function [] = display_(S,Fs,scaling_f,mute)
nbSamples = sum(size(S))-1;
top = max(S);
bottom = min(S);

plot([0:1/scaling_f:nbSamples/scaling_f-1/scaling_f], S); hold on;
ylabel('Sample value'); title('Signal');
if (scaling_f == 1)  xlabel(['Samples']); xlim([1 nbSamples]); hold on; 
    else  xlabel(['Time']); xlim([0 nbSamples/Fs]);hold on; 
end

% fs_play = 44100;

% if Fs ~= fs_play
%   playout=resample(S,fs_play,Fs);
% else
%   playout=S;
% end;
pause(0.5);
if (~mute)
    mysound(S,Fs);
%     drawnow
%     pause(0.5);
%     player = audioplayer(S, Fs);
%     play(player);
%     while player.isplaying
%         h1 = plot(get(player, 'CurrentSample')/scaling_f*ones(2), [bottom top], '-r'); hold on;
%         drawnow
%         delete(h1);
%     end
end
end

