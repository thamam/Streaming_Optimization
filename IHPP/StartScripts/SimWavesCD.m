%%  SimWavesCD
%% SimulationParamsIteration
clear ;close all; clc
addpath('BasisFunctions','NHPP_Simulation',...
    'Coord_Descent','CVX_','BwdFwdDscnt','DisplayScripts','wavesProp');
%CVX_Support_Test();
% GLOBALS - UPPER CASE VARIABLES 

global FIX_SEED
global DEBUG_MODE
% Set Global values
FIX_SEED = false;
DEBUG_MODE = false;
%% Configure simulation Settings
%Set Simulation parameters
RUNMode.CVX = true;
RUNMode.wavesCD = true;
RUNMode.CoordinateDescent = true;
lambdaMax = 25 ; %maximal rate for Lambda(t)
params = struct('Prate',lambdaMax,'splineOrder',3,'simTstart',2, 'simTfinish',59,...
    'iterThs',1e-5);
WCD_opts = struct ('sectionLength',5,'MaxIter',200);
GCD_opts = struct ('sectionLength',5,'MaxIter',500);
% loading previous sim data options0
UseSavedNhppSimdData = true;
OverWrite_NhppSimData = true;
Save_data = true;
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
%% GLOBAL SOLUTION CVX
%Estimate Entire Sample using negative log likelihood minimization
if RUNMode.CVX
    try
        if(DEBUG_MODE)
            save('DebugData.mat','-append')
        end
        run('CVX_Global_Solution_SCRIPT');
    catch M
        Exception.CVX = M;
        fprintf('CVX Exception')
    end
    params.Lambda_CVX=Lambda_CVX;
end
%% WCD
%to print results use <plotScript_coordniateDescent.m>,
% save DebugData.mat
if RUNMode.wavesCD
%     try
        clc;
        params.opts= WCD_opts;
        WCD_tic = tic;
        [WavesCoordDescOut]  = wavesCD(tK,params,WCD_opts);
        ElapsedTime_CD= toc(WCD_tic);
%     catch M
        Exception.WCD = M;
        fprintf('WCD Exception\n')
    end
end
%

%% GCD
%to print results use <plotScript_coordniateDescent.m>

if RUNMode.CoordinateDescent
    try
        clc;
        CD_tic = tic;
        [CoordDescOut]  = NHP_CoordDescent(tK,params,GCD_opts);
        ElapsedTime_CD= toc(CD_tic);
    catch M
        Exception.GCD = M;
        fprintf('GCD Exception\')
    end
end


saveRes = 'y'; %input('Would you like to save the results? (y/n)','s');
if strcmpi(saveRes,'y')
    formatOut = 'D-mm_dd_yy_T-HH_MM_PM';
    timeStamp = datestr(datetime('now'),formatOut);
    save(['results_',timeStamp,'_',params.fileNhppName]);
end


saveVideo = true;

% WCD_display_results




