function [ idx ] = sq_enc(in, n_bits, xmax, m, en_plots)
%function [ idx ] = sq_enc(in, n_bits, xmax, m)
%   
inlen = length(in);
idx = zeros(inlen,1);
L = 2^n_bits; % Number of quantization interval
delta = 2*xmax/L;%Step size
val = m-xmax+delta/2:delta:m+xmax-delta/2; %The output possible values

if en_plots
    figure, plot(in); hold on;
    for i=1:length(val)
        plot([0 inlen],val(i)*[1 1],'-r');
        hold on;
    end
    title(['USQ: Output levels. ' num2str(n_bits) ' bits, m= ' num2str(m) ', xmax= ' num2str(xmax) ]);
    xlabel('Time (Samples)'); ylabel('Signal value');
end

%To get the index we will substract one sample by all the quantization
%value and find the minimum.
columnIn = in*ones(1,length(val));%Dim inlen x L : contains  L similar columns
rowVal = ones(length(in),1)*val;%Dim inlen x L   : contains inlen similar rows 
Sub = abs(columnIn-rowVal);
[~, idx] = min(Sub,[],2); %Find the index of the minimum of each row, the result is a column.


end

