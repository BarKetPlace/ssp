function [ idx ] = sq_enc(in, n_bits, xmax, m, en_plots)
%function [ idx ] = sq_enc(in, n_bits, xmax, m)
%   
inlen = length(in);
idx = zeros(inlen,1);
L = 2^n_bits; % Number of quantization interval
delta = 2*xmax/L;
val = m-xmax:delta:m+xmax;

if en_plots
    figure, plot(in); hold on;
    for i=1:length(val)
        plot([0 inlen],val(i)*[1 1],'-r');
        hold on;
    end
    title('Quantization levels');
end


columnIn = in*ones(1,length(val));
rowVal = ones(length(in),1)*val;
Sub = abs(columnIn-rowVal);
[~, idx] = min(Sub,[],2);


% idx = interp1(in,in,val);

end

