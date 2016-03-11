function [codeA] = encodefilter(A, cb1, cb2 )
%ENCODEFILTER encode the filter parameters found in A with a 2 stages
%vector quantizer. It uses codebook 1 (cb1) to code the coefficients found
%in A and codebook 2 (cb2) to code the residual

    % get order of LP analysis
M = size(cb1) ;
M = M(2) ;

    % create codeA vector
codeA = zeros(length(A), 2) ;

    % for each polynomial vectors
for i = 1:length(A)
   
        % get lsf
    lsf = poly2lsf(A(i,:))' ;
    
        % get distance between lsf and codebook
    dist = (1/M) * sum((ones(2^10, 1)*lsf - cb1).^2, 2) ;
    
        % get index of smallest distance
    indx1 = find( dist == min(dist)) ;
    
        % get residual and code it same way as before but with lsfCB2
    res = lsf - cb1(indx1, : ) ;
    dist = (1/M) * sum((ones(2^10, 1)*res - cb2).^2, 2) ;
    indx2 = find( dist == min(dist)) ;
    
        % update codeA vector with indexes from both codebooks
    codeA(i, :) = [indx1, indx2] ;        
end
end

