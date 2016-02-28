function [outq] = sq_dec(idx, n_bits, xmax, m)
%function [outq] = sq_dec(idx, n_bits, xmax, m)
% INPUT idx : row or columns vector, 
%       n_bits : Number of bits to code de interval
%       xmax, m : the interval will be m-xmax:m+xmax
% OUTPUT outq : Row containing the quantized signal

L = 2^n_bits; % Number of quantization interval
delta = 2*xmax/L;%Step size
val = m-xmax:delta:m+xmax;%Row containing the different values

outq = val(idx); %Quantized signal



end

