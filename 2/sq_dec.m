function [outq] = sq_dec(idx, n_bits, xmax, m)
%function [outq] = sq_dec(idx, n_bits, xmax, m)
L = 2^n_bits; % Number of quantization interval
delta = 2*xmax/L;%Step size
val = m-xmax:delta:m+xmax;%Row containing the different values

outq = val(idx); %Quantized signal



end

