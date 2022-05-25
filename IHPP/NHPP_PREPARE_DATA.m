%% NHPP_PREPARE_DATA.m
%% Initialize
[params  ] = NHPP_GCD_Initialize(params); %initializing function
fs = 1e4;
tt =  params.simTstart:1/fs:params.simTfinish;
params.tt=tt;
ttEntireSupport = params.tsup.Start:1/fs:params.tsup.End;
SINTHESIZE_SIM_DATA     % Synthesize NHPP simulation Data
%% Generate random sample point of L(t)
%         MLE.theoretical = likelyHoodEstimateSpline(tK,params,alphaSyn);
 [Lambda_GroundTruth]      = lambda_t(tt,alphaSyn,params.phi_n,params.nInd);
[Lambda_at_arrivalPoints]   = lambda_t(tK,alphaSyn,params.phi_n,params.nInd);
 if(DEBUG_MODE)
    params.Lambda_GroundTruth = Lambda_GroundTruth;
    params.Lambda_at_arrivalPoints = Lambda_at_arrivalPoints ;
    params.tK=tK;
 end
if DEBUG_MODE    % Save all WS Data before running various methods
    save('DebugData.mat')
end

