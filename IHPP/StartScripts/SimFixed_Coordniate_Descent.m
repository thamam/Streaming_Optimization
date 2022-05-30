%% Sim Start - Reconstructing Lambda(t) from sample of known arrival times
%% Notation
% Lambda(t) - Ground truth intensity function
% CVX - sloves min log likelihood s.t. lambda(t)>=0
%% GLOBALS - UPPER CASE VARIABLES
global FIX_SEED
global DEBUG_MODE
% Set Global values
FIX_SEED = false;
DEBUG_MODE = true;
%% SimulationParamsIteration
clear ;close all; clc
addpath('BasisFunctions','NHPP_Simulation',...
    'Coord_Descent','CVX_','BwdFwdDscnt','DisplayScripts',...
    'wavesProp');
%CVX_Support_Test();
%% Simulation Notes
SimNotes = ' Testing end to end run of fixed GCD';
%% Configure simulation Settings
%Set Simulation parameters
RUNMode.CVX = true;
RUNMode.CoordinateDescent = true;
lambdaMax = 45 ; %maximal rate for Lambda(t)
params = struct('Prate',lambdaMax,'splineOrder',3,'simTstart',2, 'simTfinish',59,...
    'iterThs',1e-6);
GCD_opts = struct ('sectionLength',5,'MaxIter',50);
%% Multi params sim run
multiSim.LamdaMax = lambdaMax; %[25,45,55];
% multiSim.LamdaMax = [25,45,55];
for lambdaMax = multiSim.LamdaMax
    params.Prate = lambdaMax;
    % loading previous sim data options0
    UseSavedNhppSimdData = true;
    OverWrite_NhppSimData = true;
    Save_data = true;
    %% Initialize
    [params  ] = NHPP_GCD_Initialize(params); %initializing function
    fs = 1e2;
    tt =  params.simTstart:1/fs:params.simTfinish;
    ttEntireSupport = params.tsup.Start:1/fs:params.tsup.End;
    
    %Try to use simulation samples from previous runs
    if UseSavedNhppSimdData
        try
            load(params.fileNhppName);
        catch ME
            if( strcmp(ME.identifier,'MATLAB:load:couldNotReadFile') )
                ME.message
                fprintf('Simulating new process \n');
                [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
                [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
                dateOfFile=date;
                save(params.fileNhppName,'tK','alphaSyn','dateOfFile');
            else
                load(params.fileNhppName);
            end
        end
    else % Required new simulation
        fprintf('Simulating new process data \n');
        [alphaSyn] = syntContIntNHPP(params); %generate coeffs for Lambda(t)
        [ tK ] = nhppSynt(@(t)lambda_t(t,alphaSyn,params.phi_n,params.nInd) ,params); %send function handle of lambda as function of t
        if(OverWrite_NhppSimData)
            dateOfFile=date;
            save(params.fileNhppName,'tK','alphaSyn','dateOfFile');
        end
    end
    
    %% Generate random sample point of L(t)
    %         MLE.theoretical = likelyHoodEstimateSpline(tK,params,alphaSyn);
    [Lambda_GroundTruth]      = lambda_t(tt,alphaSyn,params.phi_n,params.nInd);
    [Lambda_at_arrivalPoints]   = lambda_t(tK,alphaSyn,params.phi_n,params.nInd);
    
    %% Plot theoretical lambda and it's sampled arrival times
    if false
        figure,plot(tt,Lambda_GroundTruth,'--g');title('\lambda(t)');xlabel('t'); hold all;
        scatter(tK,Lambda_at_arrivalPoints,'*k'); ylabel(' \lambda^\^(t)');xlabel('t'); xlabel('time');title(' Estimate of \lambda(t)'); legendVec={'\lambda(t)';'Events Samples'};
    end
    
    %% GLOBAL SOLUTION CVX
    %Estimate Entire Sample using negative log likelihood minimization
    if RUNMode.CVX
        run CVX_Global_Solution_SCRIPT
    end
    %% GCD
    %to print results use <plotScript_coordniateDescent.m>
    if DEBUG_MODE
        save DebugData.mat
    end
    if RUNMode.CoordinateDescent
        clc;
        CD_tic = tic;
        [CoordDescOut]  = NHP_CoordDescent(tK,params,GCD_opts);
        ElapsedTime_CD= toc(CD_tic);
    end
    %%
    saveRes = 'y'; %input('Would you like to save the results? (y/n)','s');
    if strcmpi(saveRes,'y')
        formatOut = 'D-mm_dd_yy_T-HH_MM_PM';
        timeStamp = datestr(datetime('now'),formatOut);
        save(['results_',timeStamp,'_',params.fileNhppName]);
    end
end


