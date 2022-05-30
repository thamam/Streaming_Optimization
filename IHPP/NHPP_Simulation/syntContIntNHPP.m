function [alpha] = syntContIntNHPP(params,t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    hfunc = params.phi_n;
   
   if(nargin<=2) %create time vector if required
        t=linspace(params.simTstart,params.simTfinish,100);
    end
    
    %Evaluate PSI at uniform grid 
    nInd=params.nInd; 
    PSI = zeros( numel(nInd),numel(t));
    for i=1:numel(nInd)
        PSI(i,:) = hfunc(t,nInd(i));
    end
    
    alpha = abs(randn(numel(nInd),1)); 
    alpha =  alpha*params.Prate/max(PSI.'*alpha); %scale x to achieve maximal rate 
    alpha = alpha +params.minRate;

end

