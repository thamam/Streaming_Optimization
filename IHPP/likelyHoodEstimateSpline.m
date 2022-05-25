function [ MLE ] = likelyHoodEstimateSpline(tN,params,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    phi_n = params.phi_n;
    L = params.splineOrder;
    
    %Evaluate PSI at uniform grid 
    nInd=params.nInd; 
    
    tsup = struct('Start',[],'End',[]);
    tsup.Start = params.simTstart;%-(L+1)/2;
    tsup.End = params.simTfinish;%Tsim+(L+1)/2;

    d= numel(nInd);
    
    K = numel(tN);
    c = zeros(d,1);
    for i=1:d % 
%         c(i) = integral(@(t)psi_n(t,nInd(i)),0,params.Tsim);
        c(i) = integral(@(t)phi_n(t,nInd(i)),tsup.Start,tsup.End);

    end
    
    %Assemble PSI | PSI[i,j] = psi_n(tj)
    PSI = zeros( d,K);
    for i=1:numel(nInd)
        PSI(i,:) = phi_n(tN,nInd(i));
    end
    
    MLE = alpha.'*c -ones(K,1).'*(log(PSI.'*alpha));

end

