function [  ] = nhpp_Initialize()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global params
    
    
%% General Settings
    params.fileNhppName=['NHPP_data_maxRate_',num2str(params.Prate),'_Time_',num2str(params.simTstart),...
    '_',num2str(params.simTfinish),'_SplineOrder_',num2str(params.splineOrder),'.mat'];

    params.duration = params.simTfinish - params.simTstart;
    if (params.duration <= 0 )
        error(' sim start time - sim finish time must be greater than zero');
    end
%% Randomization
    rng('default');
    rng(1);  
    params.seed = rng;
        
%% DCM Configuration Settings
    params.dcmMaxIter = 5;
    params.dcmThs=1e-4;
    params.regWeight = 0.25; %  weight of proxy function in projection minimization
%     params.beta = 0.5; 
      
%% Basis Functions  - NHPP intensity function
    bspline_L = eval(['@bspline',num2str(params.splineOrder)]);  %pick bspline* where *={1,2,3,4}
    phi_n =@(t,n) bspline_L(t-n);

    params.phi_n = phi_n; %save handle to params

%% Evaluate phi integrals vector PHI_int
    L = params.splineOrder;
    tsup.Start = params.simTstart -(L+1)/2;
    tsup.End = params.simTfinish + (L+1)/2;
    params.tsup=tsup;
%     nInd=ceil(tsup.Start):floor(tsup.End); 
    [ nInd ] = findSupportIndexes( phi_n,tsup,params.simTstart,params.simTfinish );
    params.nInd = nInd;
    d= numel(nInd);   
    PHI_int = zeros(d,1);
    for i=1:d % 
        PHI_int(i) = integral(@(t)phi_n(t,nInd(i)),tsup.Start,tsup.End);
    end        
    
%% Assemble PHI
%PHI:= [phi_n(t-0),phi_n(t-1),...,phi_n(t-N)].' function handle as function of t
    PHI = @(t) phi_n(t,nInd).';
%% Assemble negative-log-likelihood's gradient handle function
    params.h_g = @(alpha,tK) PHI_int - sum(logsumPHI (alpha,tK),2);

    function [logPHImat] = logsumPHI(alpha,tK)
        logPHImat = zeros(numel(alpha),numel(tK));
        for k=1:numel(tK)
            logPHImat(:,k) = (1/(alpha.'*PHI(tK(k))))*PHI(tK(k));
        end       
    end

end

