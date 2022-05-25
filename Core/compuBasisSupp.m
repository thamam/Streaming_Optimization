function [minTsup, maxTsup] = compuBasisSupp(tsect, spord, nInd)
    % the function computes the maximal time support for basis indexed by
    % nInd with hard limit given by tsect.
L = (spord+1 )/ 2; %One side of the support
    
    minTsup = max(min(nInd) - L, tsect(1));
    maxTsup = min(max(nInd) + L, tsect(2));
%     maxbasissupp = [minTsup, maxTsup ] ;
    
%     nindSuppPrm=[];
%     nindSuppTail=[];
%     support = tsect + [-L, +L] ;
%     nindSupp = nInd( nInd>support(1) & nInd<support(2)) ;
%     nindSuppPrm = nindSupp(:) ;

end