clear all
close all
clc


load assignment1.mat

S = male_short;
Fs = 8000;
scaling_f = 1;
mute = 1;

x = S;
%Frame length 
flen=30;%ms
alen = flen/1000*Fs;
voicedthreshold = .4; %Thresholdvoiced/unvoiced
naf = ceil(length(x)/alen);

ulen =1;
M = 9;

[E, ZC, V, A,P]=analysis(x,alen,ulen,voicedthreshold,M);
%% 
close all
figure, plot(ZC); hold on; stem(V*max(ZC(:)));


figure(1);clf;
subplot(3,2,1)
plot(x) % Plot the input waveform
axis([1 length(x) min(x) max(x)]);
title('Waveform');

subplot(3,2,2)
plot(sqrt(E)) % Plot the standard deviation
axis([1 length(E) min(sqrt(E)) max(sqrt(E))]);
title('Standard deviation');
subplot(3,2,3)
plot(V) % Plot voiced/unvoiced decision
axis([1 length(V) 0 1]);
title('Voiced (1) unvoiced(0)');

subplot(3,2,4)
plot(ZC); hold on;% Plot the normalized number of zero-crossings
plot(length(ZC)*[0 1], voicedthreshold*[1 1], '-r');

axis([1 length(ZC) min(ZC) max(ZC)]);
title('Normalized number of ZC');

subplot(3,2,5)
F = Fs./P;
plot(F) % Plot the fundamental frequency in Hz
% axis([1 length(F) 0 600]); 
ylabel('Frequency (Hz)');
title('Fundamental frequency');

subplot(3,2,6)
S = zeros(512, naf);
for n=1:naf
S(:,n) = 20*log10(abs(freqz(1,A(n,:),512)));
end
S = flipud(S);
colormap(gray);
imagesc(S); % Illustrate the vocal tract envelope in a spectrogram
%style!
title('Vocal tract envelope');

