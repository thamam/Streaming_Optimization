function [params  ] = NHPP_GCD_Initialize(params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here    
global FIX_SEED    
%% General Settings
    params.fileNhppName=['GCD_simData_maxRate_',num2str(params.Prate),'_T_',num2str(params.simTstart),...
    '_',num2str(params.simTfinish),'_L_',num2str(params.splineOrder),'.mat'];

    params.duration = params.simTfinish - params.simTstart;
    if (params.duration <= 0 )
        error(' sim start time - sim finish time must be greater than zero');
    end
%% Randomization
if (FIX_SEED)
    rng('default');
    rng(1);  
    params.seed = rng;
end
%% Basis Functions  - NHPP intensity function
    bspline_L = eval(['@bspline',num2str(params.splineOrder)]);  %pick bspline* where *={1,2,3,4}
    phi_n =@(t,n) bspline_L(t-n);
    params.phi_n = phi_n; %save handle to params

%% Evaluate phi integrals vector PHI_int
    L = params.splineOrder;
    [ nInd ] = findSupportIndexes( phi_n,params.simTstart,params.simTfinish,L );
    params.SupLen = (L+1)/2;
    tsup.Start = max(params.simTstart, nInd(1) -(L+1)/2);
    tsup.End = min( params.simTfinish, nInd(end) + (L+1)/2);
    params.tsup=tsup;
    params.nInd = nInd;
    d= numel(nInd);   
    params.d=d;
%     PHI_int = zeros(d,1);
%     for i=1:d % 
%         PHI_int(i) = integral(@(t)phi_n(t,nInd(i)),tsup.Start,tsup.End);
%     end        
    
%% Assemble PHI
%PHI:= [phi_n(t-0),phi_n(t-1),...,phi_n(t-N)].' function handle as function of t
%     PHI = @(t) phi_n(t,nInd).';

end

