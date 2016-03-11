function [Aq] = decodefilter(codeA, cb1, cb2)
%DECODEFILTER decode into LP parameters coefficients the indexes found in
%codeA using codebooks cb1 and cb2.

    % create vector Aq
M = size(cb1) ;
M = M(2) + 1 ;
Aq = zeros(length(codeA), M) ;

for i = 1:length(codeA)
   
        % get indexes
    indxs = codeA(i, :) ;
    
        % get lsf et residual
    lsf = cb1(indxs(1), :) ;
    res = cb2(indxs(2), :) ;
    lsf = sort(lsf + res) ;         % use sort to make sure it represents a minimum phase whitening filter

        % call lsf2poly and store the result
    Aq(i, :) = lsf2poly(lsf) ;
end
end

