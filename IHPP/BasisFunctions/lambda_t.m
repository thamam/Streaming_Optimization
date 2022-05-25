function [Lambda ] = lambda_t(t,alpha_hat,phi_n,nInd )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here  

    d= numel(nInd);
    K=numel(t);
    PSI = zeros( d,K);
    for i=1:numel(nInd)
        PSI(i,:) = phi_n(t,nInd(i));
    end
        
    Lambda = PSI.'*alpha_hat;
    
    if (size(t,1)~= size(Lambda,1))
        Lambda = Lambda.';
    end
end

