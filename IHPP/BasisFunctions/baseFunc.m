function [ f ] = baseFunc( t,alpha,J )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    recFunc = @(t) 1*((t>=0) & (t<1));
    recFunc2j = @(t,n,j) 2^(j /2)* recFunc(2^j*t-n);
    %%
    d = numel(alpha);
    nVec = 0:d-1;
    Xtemp = zeros(d,numel(t));
    for q=1:numel(nVec)
        n=nVec(q);
        Xtemp(q,:) = alpha(q)*recFunc2j(t,n,J);      
    end
       
    f = sum(Xtemp,1);

end

